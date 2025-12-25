// src/Sparrow.cpp
#include "Sparrow.h"
#include "Collision.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <chrono>
#include <random>
#include <iomanip>

SeparatorWorker::SeparatorWorker(int _id, int seed, const std::vector<CompositeShape>& initItems, size_t n) 
    : id(_id), items(initItems), rng(seed) {
    weights.resize(n, std::vector<double>(n, 1.0));
    current_energy = 0.0;
}

void SeparatorWorker::load(const std::vector<CompositeShape>& masterItems, const std::vector<std::vector<double>>& masterWeights) {
    // Tái sử dụng bộ nhớ vector, chỉ copy dữ liệu
    items = masterItems;
    weights = masterWeights;
    // Không cần tính lại energy ngay, sẽ tính sau khi move
}

double SeparatorWorker::evaluateSample(int itemIdx, Vec2 pos, double angle, const SparrowConfig& config) {
    CompositeShape tempItem = items[itemIdx];
    tempItem.setTransform(pos, angle); 
    
    double totalCost = 0.0;
    double halfSize = config.container_size / 2.0;

    // Wall Penalty
    const AABB& aabb = tempItem.totalAABB;
    double outOfBounds = 0.0;
    if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
    if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
    if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
    if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
    
    if (outOfBounds > 0) {
        totalCost += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
    }

    // Overlap
    for (size_t j = 0; j < items.size(); ++j) {
        if ((int)j == itemIdx) continue;
        if (tempItem.totalAABB.overlaps(items[j].totalAABB)) {
             double overlap = quantify_collision(tempItem, items[j]);
             if (overlap > 1e-9) {
                 totalCost += overlap * weights[itemIdx][j];
             }
        }
    }
    return totalCost;
}

void SeparatorWorker::searchPosition(int itemIdx, const SparrowConfig& config) {
    const CompositeShape& current = items[itemIdx];
    double current_e = evaluateSample(itemIdx, current.pos, current.angle, config);

    // Static thread_local để tái sử dụng memory cho candidates trong luồng này
    static thread_local std::vector<std::pair<double, Transform>> candidates;
    candidates.clear();
    candidates.reserve(config.n_samples);

    double halfSize = config.container_size / 2.0;
    double foc_radius = config.container_size * 0.15;
    
    std::uniform_real_distribution<double> dist01(0.0, 1.0);
    std::uniform_real_distribution<double> distContainer(-halfSize, halfSize);
    std::uniform_real_distribution<double> distFoc(-foc_radius, foc_radius);
    std::uniform_real_distribution<double> distAngle(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> distJitter(-0.02, 0.02);

    int n_div = (int)(config.n_samples * 0.4);
    int n_foc = (int)(config.n_samples * 0.4);
    int n_flip = config.n_samples - n_div - n_foc;

    // Phase 1: Global
    for (int s = 0; s < n_div; ++s) {
        double x = distContainer(rng);
        double y = distContainer(rng);
        double rot = (dist01(rng) < 0.5) ? 0.0 : PI; 
        rot += distJitter(rng);
        candidates.emplace_back(evaluateSample(itemIdx, {x, y}, rot, config), Transform({x, y}, rot));
    }

    // Phase 2: Local
    for (int s = 0; s < n_foc; ++s) {
        double x = current.pos.x + distFoc(rng);
        double y = current.pos.y + distFoc(rng);
        double rot = current.angle + distJitter(rng);
        candidates.emplace_back(evaluateSample(itemIdx, {x, y}, rot, config), Transform({x, y}, rot));
    }

    // Phase 3: Flip
    for (int s = 0; s < n_flip; ++s) {
        double x = current.pos.x + distFoc(rng);
        double y = current.pos.y + distFoc(rng);
        double rot = current.angle + PI + distJitter(rng);
        while (rot > 2*PI) rot -= 2*PI;
        candidates.emplace_back(evaluateSample(itemIdx, {x, y}, rot, config), Transform({x, y}, rot));
    }

    // Sort & Pick Best
    std::sort(candidates.begin(), candidates.end());

    // Coordinate Descent (Local Refinement)
    int K = 3; 
    Transform best_t(current.pos, current.angle);
    double best_e = current_e;

    for (int i = 0; i < std::min(K, (int)candidates.size()); ++i) {
        Transform t = candidates[i].second;
        double e = candidates[i].first;
        if (e > current_e * 1.5) continue;

        double step = foc_radius * 0.5;
        double rot_step = 0.05; 
        
        for (int iter = 0; iter < 10; ++iter) { 
             bool improved = false;
             for (double d : {-step, step}) {
                 double ne = evaluateSample(itemIdx, t.pos + Vec2{d,0}, t.angle, config);
                 if (ne < e) { t.pos.x += d; e = ne; improved = true; }
             }
             for (double d : {-step, step}) {
                 double ne = evaluateSample(itemIdx, t.pos + Vec2{0,d}, t.angle, config);
                 if (ne < e) { t.pos.y += d; e = ne; improved = true; }
             }
             for (double d : {-rot_step, rot_step}) {
                 double ne = evaluateSample(itemIdx, t.pos, t.angle + d, config);
                 if (ne < e) { t.angle += d; e = ne; improved = true; }
             }
             if (!improved) { step *= 0.5; rot_step *= 0.5; }
             if (step < 1e-4) break;
        }
        if (e < best_e) { best_e = e; best_t = t; }
    }

    if (best_e < current_e) {
        items[itemIdx].setTransform(best_t.pos, best_t.angle);
    }
}

double SeparatorWorker::calculate_total_energy(const SparrowConfig& config) {
    double energy = 0.0;
    double halfSize = config.container_size / 2.0;

    for (size_t i = 0; i < items.size(); ++i) {
        const AABB& aabb = items[i].totalAABB;
        double outOfBounds = 0.0;
        if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
        if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
        if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
        if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
        if (outOfBounds > 0) {
            energy += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
        }

        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[i].totalAABB.overlaps(items[j].totalAABB)) {
                energy += quantify_collision(items[i], items[j]);
            }
        }
    }
    return energy;
}

void SeparatorWorker::move_items(const SparrowConfig& config) {
    // 1. Identify Colliding Items
    std::vector<int> Ic;
    double halfSize = config.container_size / 2.0;
    
    for (int i = 0; i < (int)items.size(); ++i) {
        bool colliding = false;
        const AABB& aabb = items[i].totalAABB;
        if (aabb.min.x < -halfSize || aabb.max.x > halfSize || 
            aabb.min.y < -halfSize || aabb.max.y > halfSize) {
            colliding = true;
        } else {
            for (int j = 0; j < (int)items.size(); ++j) {
                if (i == j) continue;
                if (items[i].totalAABB.overlaps(items[j].totalAABB)) {
                    if (quantify_collision(items[i], items[j]) > 1e-9) {
                        colliding = true; break;
                    }
                }
            }
        }
        if (colliding) Ic.push_back(i);
    }

    if (Ic.empty()) {
        current_energy = 0.0;
        return;
    }

    // Shuffle items to move
    std::shuffle(Ic.begin(), Ic.end(), rng);

    // Try to move colliding items
    for (int idx : Ic) {
        searchPosition(idx, config);
    }

    // Recalculate energy after moves
    current_energy = calculate_total_energy(config);
}


// ==========================================
// SPARROW SOLVER IMPLEMENTATION
// ==========================================

// Helper
double SparrowSolver::randomDouble(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

void SparrowSolver::saveState() {
    savedItems = items;
    savedContainerSize = config.container_size;
}

void SparrowSolver::loadState() {
    items = savedItems;
    config.container_size = savedContainerSize;
}

// Hàm khởi tạo cũ (Giữ nguyên logic)
void SparrowSolver::initializePlacement() {
    double halfSize = config.container_size / 2.0;
    int n = items.size();
    int rows = std::ceil(std::sqrt(n));
    double spacing = config.container_size / rows;
    int idx = 0;
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < rows && idx < n; ++c) {
            int itemIdx = indices[idx];
            double x = -halfSize + spacing * (c + 0.5);
            double y = -halfSize + spacing * (r + 0.5);
            double rot = ((r + c) % 2 == 0) ? 0.0 : PI; 
            rot += randomDouble(-0.1, 0.1); 
            double jitter = spacing * 0.1;
            x += randomDouble(-jitter, jitter);
            y += randomDouble(-jitter, jitter);
            items[itemIdx].setTransform({x, y}, rot);
            idx++;
        }
    }
}

SparrowSolver::SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg) 
    : items(inputItems), config(cfg) {

    // 1. Config Threads
    if (config.n_threads <= 0) {
        config.n_threads = omp_get_num_procs();
    }
    omp_set_num_threads(config.n_threads);

    std::random_device rd;
    rng = std::mt19937(rd());

    // 2. Init Weights
    size_t n = items.size();
    weights.resize(n, std::vector<double>(n, 1.0));
    for (auto& item : items) item.precomputeAllPoles(0.01);
    
    initializePlacement();

    // 3. CREATE PERSISTENT WORKERS
    workers.reserve(config.n_threads);
    for (int i = 0; i < config.n_threads; ++i) {
        // Mỗi worker có seed riêng để chạy độc lập
        workers.emplace_back(i, rd() + i, items, n);
    }
    
    std::cout << ">>> Solver Initialized with " << config.n_threads << " persistent workers." << std::endl;
}

double SparrowSolver::getTotalEnergy() const {
    // Tính toán trên Master
    double energy = 0.0;
    double halfSize = config.container_size / 2.0;
    for (size_t i = 0; i < items.size(); ++i) {
        const AABB& aabb = items[i].totalAABB;
        double outOfBounds = 0.0;
        if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
        if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
        if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
        if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
        if (outOfBounds > 0) energy += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[i].totalAABB.overlaps(items[j].totalAABB)) {
                energy += quantify_collision(items[i], items[j]);
            }
        }
    }
    return energy;
}

void SparrowSolver::updateMasterWeights() {
    double e_max = 0.0;
    size_t n = items.size();
    for (size_t a = 0; a < n; ++a) {
        for (size_t b = a + 1; b < n; ++b) {
            if (items[a].totalAABB.overlaps(items[b].totalAABB)) {
                double e = quantify_collision(items[a], items[b]);
                if (e > e_max) e_max = e;
            }
        }
    }
    
    if (e_max < 1e-9) return;

    for (size_t a = 0; a < n; ++a) {
        for (size_t b = a + 1; b < n; ++b) {
            double e = 0.0;
            if (items[a].totalAABB.overlaps(items[b].totalAABB)) {
                e = quantify_collision(items[a], items[b]);
            }
            double m = (e > 1e-9) ? (config.Ml + (config.Mu - config.Ml) * (e / e_max)) : config.Md;
            weights[a][b] = std::max(1.0, weights[a][b] * m);
            weights[b][a] = weights[a][b]; 
        }
    }
}

// === ALGORITHM 9: SEPARATE (REFACTORED LIKE RUST) ===
bool SparrowSolver::separate(int kmax, int nmax, std::chrono::steady_clock::time_point deadline) {
    double best_global_energy = getTotalEnergy();
    if (best_global_energy <= 1e-5) return true;
    if (std::chrono::steady_clock::now() > deadline) return false;

    int k = 0;
    while (k < kmax) {
        int n = 0;
        bool improved_in_try = false;

        while (n < nmax) {
            if (std::chrono::steady_clock::now() > deadline) return false;

            // 1. SYNC & PARALLEL EXECUTION (Fork-Join)
            // Không có critical section bên trong vòng lặp của Worker!
            #pragma omp parallel for
            for (int i = 0; i < (int)workers.size(); ++i) {
                // Load Master State
                workers[i].load(items, weights);
                // Worker tự chạy logic di chuyển của riêng mình
                workers[i].move_items(config);
            }

            // 2. GATHER RESULT (Master Thread)
            int best_worker_idx = -1;
            double min_worker_energy = std::numeric_limits<double>::max();

            for (int i = 0; i < (int)workers.size(); ++i) {
                if (workers[i].current_energy < min_worker_energy) {
                    min_worker_energy = workers[i].current_energy;
                    best_worker_idx = i;
                }
            }

            // 3. UPDATE MASTER STATE
            if (best_worker_idx != -1) {
                // Cập nhật Master theo worker giỏi nhất
                items = workers[best_worker_idx].items;
                // Lưu ý: Rust update weights ở bước sau, không copy weights từ worker về
                
                if (min_worker_energy < best_global_energy) {
                    best_global_energy = min_worker_energy;
                    improved_in_try = true;
                }
            }

            if (best_global_energy <= 1e-5) return true;

            // 4. UPDATE WEIGHTS (Algorithm 8)
            // Chỉ Master làm việc này 1 lần sau mỗi bước sync
            updateMasterWeights();

            // Stagnation logic
            if (min_worker_energy < best_global_energy * 0.999) { // Có cải thiện đáng kể
                n = 0;
            } else {
                n++;
            }
        }
        
        k++;
        if (improved_in_try) k = 0; // Nếu strike này có tiến bộ, reset strike
    }

    return (best_global_energy <= 1e-5);
}

// Scramble logic cũng cần cập nhật
void SparrowSolver::scrambleItems() {
    // ... Logic giống cũ nhưng thực hiện trên Master items ...
    int count = std::max(1, (int)(items.size() * 0.2));
    std::vector<int> indices(items.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    double halfSize = config.container_size / 2.0;

    for (int i = 0; i < count; ++i) {
        int idx = indices[i];
        items[idx].setTransform(
            {randomDouble(-halfSize, halfSize), randomDouble(-halfSize, halfSize)},
            (randomDouble(0,1) < 0.5) ? (items[idx].angle + PI) : randomDouble(0, 2*PI)
        );
    }
}

// Các hàm explore, compress, descent, solve giữ nguyên logic điều khiển luồng (Flow Control)
// nhưng gọi separate() mới.
// ... (Copy lại logic explore/compress/solve từ code trước của bạn) ...
void SparrowSolver::explore() {
    auto start_time = std::chrono::steady_clock::now();
    auto deadline = start_time + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(config.TLx));
    
    if (!separate(config.Kx, config.Nx * 2, deadline)) {
        if (std::chrono::steady_clock::now() > deadline) return;
        config.container_size *= 1.05;
        initializePlacement(); 
        separate(config.Kx, config.Nx * 2, deadline);
    }

    double best_size = config.container_size;
    std::vector<CompositeShape> s_best = items;

    int iter = 0;
    int fail_streak = 0;
    double current_Rx = config.Rx; 

    while (iter < config.max_outer_loops) {
        if (std::chrono::steady_clock::now() > deadline) break;
        
        if (current_Rx < 1e-5) { current_Rx = config.Rx; scrambleItems(); }

        std::vector<CompositeShape> s_prev_feasible = items;
        double size_prev = config.container_size;

        config.container_size *= (1.0 - current_Rx);
        std::cout << "Explore [" << iter << "] | Rx: " << current_Rx << " | Target: " << config.container_size << std::endl;

        double scale = config.container_size / size_prev;
        for(auto& it : items) it.setTransform(it.pos * scale, it.angle);
        for(auto& r : weights) std::fill(r.begin(), r.end(), 1.0); // Reset weights

        bool feasible = separate(config.Kx, config.Nx + (fail_streak*50), deadline);

        if (feasible) {
            best_size = config.container_size;
            s_best = items;
            fail_streak = 0;
        } else {
            if (std::chrono::steady_clock::now() > deadline) {
                items = s_prev_feasible; config.container_size = size_prev; break;
            }
            fail_streak++;
            items = s_prev_feasible;
            config.container_size = size_prev;
            current_Rx *= 0.5;
            scrambleItems(); 
            std::cout << "   -> Fail. Revert & Reduce Rx to " << current_Rx << std::endl;
        }
        iter++;
    }
    items = s_best;
    config.container_size = best_size;
}

void SparrowSolver::compress(double ratio) {
    // Logic compress giữ nguyên, chỉ đảm bảo gọi separate mới
    auto start_time = std::chrono::steady_clock::now();
    auto deadline = start_time + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(config.TLc));
    std::vector<CompositeShape> s_star = items;
    double best_size = config.container_size;
    
    int loop_iter = 0;
    while (loop_iter < config.max_outer_loops) {
        if (std::chrono::steady_clock::now() > deadline) break;
        std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time;
        double progress = elapsed.count() / config.TLc;
        if (progress > 1.0) progress = 1.0;
        double r = config.Rs_c + (config.Re_c - config.Rs_c) * progress;

        items = s_star;
        double old_size = config.container_size;
        config.container_size *= (1.0 - r);
        for(auto& row : weights) std::fill(row.begin(), row.end(), 1.0);

        if (separate(config.Kc, config.Nc, deadline)) {
            s_star = items; best_size = config.container_size;
            std::cout << "   Compress success -> " << best_size << std::endl;
        } else {
            config.container_size = old_size;
            std::cout << "   Compress fail." << std::endl;
        }
        loop_iter++;
    }
    items = s_star; config.container_size = best_size;
}

void SparrowSolver::descent() { compress(0.0); }

void SparrowSolver::solve() {
    // Logic solve giữ nguyên
    std::cout << ">>> Solver Started. Items: " << items.size() << std::endl;
    int attempts = 0;
    while (attempts < 10) {
        auto deadline = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(config.TLx));
        if (separate(config.Kx, config.Nx * 2, deadline)) break;
        config.container_size *= 1.1; initializePlacement(); attempts++;
    }
    if (attempts >= 10) return;

    double bestSize = config.container_size;
    std::vector<CompositeShape> bestItems = items;
    
    for (int iter = 0; iter < config.max_restart; ++iter) {
        std::cout << "\n=== Major Iteration " << iter << " ===" << std::endl;
        descent();
        if (config.container_size < bestSize) { bestSize = config.container_size; bestItems = items; }
        items = bestItems; config.container_size = bestSize;
        explore();
        if (config.container_size < bestSize) { bestSize = config.container_size; bestItems = items; }
    }
    items = bestItems; config.container_size = bestSize;
    std::cout << ">>> Final Size: " << std::setprecision(10) << std::fixed << config.container_size << std::endl;
}