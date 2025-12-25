#include "Sparrow.h"
#include "Collision.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>

// ==========================================
// CONFIG LOADER
// ==========================================
void SparrowConfig::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return;
    std::string line;
    while (std::getline(file, line)) {
        size_t c = line.find('#');
        if (c != std::string::npos) line = line.substr(0, c);
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string k = line.substr(0, eq);
        std::string v = line.substr(eq + 1);
        auto trim = [](std::string& s) { 
            s.erase(0, s.find_first_not_of(" \t\r\n")); 
            s.erase(s.find_last_not_of(" \t\r\n") + 1); 
        };
        trim(k); trim(v);
        try {
            if (k == "container_size") container_size = std::stod(v);
            else if (k == "n_samples") n_samples = std::stoi(v);
            else if (k == "n_threads") n_threads = std::stoi(v);
            else if (k == "Nc") Nc = std::stoi(v);
            else if (k == "Nx") Nx = std::stoi(v);
            else if (k == "Kc") Kc = std::stoi(v);
            else if (k == "Kx") Kx = std::stoi(v);
            else if (k == "TLx") TLx = std::stod(v);
            else if (k == "TLc") TLc = std::stod(v);
            else if (k == "Rx") Rx = std::stod(v);
            else if (k == "Rs_c") Rs_c = std::stod(v);
            else if (k == "Re_c") Re_c = std::stod(v);
            else if (k == "max_outer_loops") max_outer_loops = std::stoi(v);
            else if (k == "max_restart") max_restart = std::stoi(v);
        } catch (...) {}
    }
}

// ==========================================
// SEPARATOR WORKER IMPLEMENTATION
// ==========================================

SeparatorWorker::SeparatorWorker(int _id, int seed, const std::vector<CompositeShape>& initItems, size_t n) 
    : id(_id), items(initItems), rng(seed) {
    // Khởi tạo bộ nhớ cục bộ
    weights.resize(n, std::vector<double>(n, 1.0));
    current_energy = 0.0;
    // Reserve trước bộ nhớ cho candidates để không bao giờ phải malloc trong loop
    candidate_buffer.reserve(4096); 
}

void SeparatorWorker::load(const std::vector<CompositeShape>& masterItems, const std::vector<std::vector<double>>& masterWeights) {
    // Tái sử dụng capacity vector, chỉ copy dữ liệu (Rất nhanh)
    items = masterItems;
    weights = masterWeights;
    // Chưa cần tính energy ngay, vì move_items sẽ tính lại
}

double SeparatorWorker::randomDouble(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

double SeparatorWorker::evaluateSample(int itemIdx, Vec2 pos, double angle, const SparrowConfig& config) {
    // Tạo bản sao tạm trên stack (rất nhẹ)
    CompositeShape tempItem = items[itemIdx];
    tempItem.setTransform(pos, angle); 
    
    double totalCost = 0.0;
    double halfSize = config.container_size / 2.0;

    // 1. Wall Penalty
    const AABB& aabb = tempItem.totalAABB;
    double outOfBounds = 0.0;
    if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
    if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
    if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
    if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
    
    if (outOfBounds > 0) {
        // Phạt lũy thừa để ép vật thể vào trong
        totalCost += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
    }

    // 2. Overlap với các vật thể khác (Sử dụng dữ liệu cục bộ của Worker)
    // Dùng Broadphase AABB trước
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

    // Xóa buffer cũ, nhưng giữ nguyên vùng nhớ đã cấp phát
    candidate_buffer.clear();
    
    double halfSize = config.container_size / 2.0;
    double foc_radius = config.container_size * 0.15; // Vùng tìm kiếm cục bộ

    int n_div = (int)(config.n_samples * 0.4); // Global search
    int n_foc = (int)(config.n_samples * 0.4); // Local search
    int n_flip = config.n_samples - n_div - n_foc; // Flip search

    // Phase 1: Global Random (T_div)
    for (int s = 0; s < n_div; ++s) {
        double x = randomDouble(-halfSize, halfSize);
        double y = randomDouble(-halfSize, halfSize);
        // Bias rotation: Ưu tiên 0 và PI cho hình đối xứng
        double rot = (randomDouble(0, 1) < 0.5) ? 0.0 : PI; 
        rot += randomDouble(-0.1, 0.1); // Jitter
        
        double e = evaluateSample(itemIdx, {x, y}, rot, config);
        candidate_buffer.emplace_back(e, Transform({x, y}, rot));
    }

    // Phase 2: Local Search (T_foc)
    for (int s = 0; s < n_foc; ++s) {
        double x = current.pos.x + randomDouble(-foc_radius, foc_radius);
        double y = current.pos.y + randomDouble(-foc_radius, foc_radius);
        double rot = current.angle + randomDouble(-0.1, 0.1);
        double e = evaluateSample(itemIdx, {x, y}, rot, config);
        candidate_buffer.emplace_back(e, Transform({x, y}, rot));
    }

    // Phase 3: Flip Strategy (Quan trọng cho vật đối xứng)
    for (int s = 0; s < n_flip; ++s) {
        double x = current.pos.x + randomDouble(-foc_radius, foc_radius);
        double y = current.pos.y + randomDouble(-foc_radius, foc_radius);
        // Lật ngược 180 độ
        double rot = current.angle + PI + randomDouble(-0.1, 0.1);
        while (rot > 2*PI) rot -= 2*PI;
        
        double e = evaluateSample(itemIdx, {x, y}, rot, config);
        candidate_buffer.emplace_back(e, Transform({x, y}, rot));
    }

    // Sắp xếp chọn ứng viên tốt nhất
    std::sort(candidate_buffer.begin(), candidate_buffer.end());

    // Local Refinement (Coordinate Descent) cho Top 3
    int K = 3; 
    Transform best_t(current.pos, current.angle);
    double best_e = current_e;

    for (int i = 0; i < std::min(K, (int)candidate_buffer.size()); ++i) {
        Transform t = candidate_buffer[i].second;
        double e = candidate_buffer[i].first;
        if (e > current_e * 1.5) continue; // Bỏ qua nếu quá tệ

        double step = foc_radius * 0.5;
        double rot_step = 0.05; 
        
        // Tinh chỉnh từng bước nhỏ
        for (int iter = 0; iter < 10; ++iter) { 
             bool improved = false;
             // X
             for (double d : {-step, step}) {
                 double ne = evaluateSample(itemIdx, t.pos + Vec2{d,0}, t.angle, config);
                 if (ne < e) { t.pos.x += d; e = ne; improved = true; }
             }
             // Y
             for (double d : {-step, step}) {
                 double ne = evaluateSample(itemIdx, t.pos + Vec2{0,d}, t.angle, config);
                 if (ne < e) { t.pos.y += d; e = ne; improved = true; }
             }
             // Angle
             for (double d : {-rot_step, rot_step}) {
                 double ne = evaluateSample(itemIdx, t.pos, t.angle + d, config);
                 if (ne < e) { t.angle += d; e = ne; improved = true; }
             }
             if (!improved) { step *= 0.5; rot_step *= 0.5; }
             if (step < 1e-5) break;
        }
        if (e < best_e) { best_e = e; best_t = t; }
    }

    // Nếu tìm được vị trí tốt hơn, cập nhật trạng thái Worker
    if (best_e < current_e) {
        items[itemIdx].setTransform(best_t.pos, best_t.angle);
    }
}

double SeparatorWorker::calculate_total_energy(const SparrowConfig& config) {
    double energy = 0.0;
    double halfSize = config.container_size / 2.0;

    for (size_t i = 0; i < items.size(); ++i) {
        // Wall Energy
        const AABB& aabb = items[i].totalAABB;
        double outOfBounds = 0.0;
        if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
        if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
        if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
        if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
        if (outOfBounds > 0) {
            energy += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
        }
        
        // Collision Energy
        for (size_t j = i + 1; j < items.size(); ++j) {
            if (items[i].totalAABB.overlaps(items[j].totalAABB)) {
                energy += quantify_collision(items[i], items[j]);
            }
        }
    }
    return energy;
}

void SeparatorWorker::move_items(const SparrowConfig& config) {
    // 1. Xác định các item đang va chạm
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

    // Shuffle thứ tự xử lý
    std::shuffle(Ic.begin(), Ic.end(), rng);

    // 2. Thử di chuyển từng item va chạm
    for (int idx : Ic) {
        searchPosition(idx, config);
    }

    // 3. Tính lại energy tổng kết của Worker này sau khi di chuyển
    current_energy = calculate_total_energy(config);
}


// ==========================================
// SPARROW SOLVER IMPLEMENTATION (MASTER)
// ==========================================

void SparrowSolver::saveState() {
    savedItems = items;
    savedContainerSize = config.container_size;
}

void SparrowSolver::loadState() {
    items = savedItems;
    config.container_size = savedContainerSize;
}

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
            rot += std::uniform_real_distribution<double>(-0.1, 0.1)(rng); 
            double jitter = spacing * 0.1;
            x += std::uniform_real_distribution<double>(-jitter, jitter)(rng);
            y += std::uniform_real_distribution<double>(-jitter, jitter)(rng);
            items[itemIdx].setTransform({x, y}, rot);
            idx++;
        }
    }
}

SparrowSolver::SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg) 
    : items(inputItems), config(cfg) {

    // 1. Cấu hình luồng
    if (config.n_threads <= 0) {
        config.n_threads = omp_get_num_procs();
    }
    omp_set_num_threads(config.n_threads);

    std::random_device rd;
    rng = std::mt19937(rd());

    // 2. Khởi tạo Weights
    size_t n = items.size();
    weights.resize(n, std::vector<double>(n, 1.0));
    for (auto& item : items) item.precomputeAllPoles(0.01);
    
    initializePlacement();

    // 3. TẠO CÁC CÔNG NHÂN KIÊN ĐỊNH (Persistent Workers)
    // Mỗi worker được gán một ID và một seed ngẫu nhiên khác nhau
    workers.reserve(config.n_threads);
    for (int i = 0; i < config.n_threads; ++i) {
        workers.emplace_back(i, rd() + i, items, n);
    }
    
    std::cout << ">>> Solver Initialized with " << config.n_threads << " persistent workers." << std::endl;
}

double SparrowSolver::getTotalEnergy() const {
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
    // Tìm max overlap trên Master state
    for (size_t a = 0; a < n; ++a) {
        for (size_t b = a + 1; b < n; ++b) {
            if (items[a].totalAABB.overlaps(items[b].totalAABB)) {
                double e = quantify_collision(items[a], items[b]);
                if (e > e_max) e_max = e;
            }
        }
    }
    
    if (e_max < 1e-9) return;

    // Cập nhật weights
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

// === ALGORITHM 9: SEPARATE (RUST STYLE) ===
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
            // - Master copy state sang Worker (trong load)
            // - Worker chạy độc lập (trong move_items)
            #pragma omp parallel for
            for (int i = 0; i < (int)workers.size(); ++i) {
                // Workers load items & weights từ Master (Deep copy nhưng an toàn)
                workers[i].load(items, weights);
                // Worker chạy logic của riêng mình
                workers[i].move_items(config);
            }

            // 2. GATHER RESULT (Master Thread)
            // Sau khi tất cả worker xong việc (barrier của parallel for), Master kiểm tra ai giỏi nhất
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
                // Master chấp nhận trạng thái của worker tốt nhất
                items = workers[best_worker_idx].items;
                
                if (min_worker_energy < best_global_energy) {
                    best_global_energy = min_worker_energy;
                    improved_in_try = true;
                    // std::cout << "   [Sep] Improved: " << best_global_energy << std::endl;
                }
            }

            if (best_global_energy <= 1e-5) return true;

            // 4. UPDATE WEIGHTS (Algorithm 8)
            // Chỉ Master tính toán lại trọng số dựa trên trạng thái mới nhất
            updateMasterWeights();

            // Stagnation logic
            if (min_worker_energy < best_global_energy * 0.999) { 
                n = 0; // Reset nếu có cải thiện đáng kể
            } else {
                n++;
            }
        }
        
        k++;
        if (improved_in_try) k = 0; 
    }

    return (best_global_energy <= 1e-5);
}

void SparrowSolver::scrambleItems() {
    int count = std::max(1, (int)(items.size() * 0.3)); // Tăng tỷ lệ scramble
    std::vector<int> indices(items.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    double halfSize = config.container_size / 2.0;

    for (int i = 0; i < count; ++i) {
        int idx = indices[i];
        // Teleport
        items[idx].setTransform(
            {std::uniform_real_distribution<double>(-halfSize, halfSize)(rng), 
             std::uniform_real_distribution<double>(-halfSize, halfSize)(rng)},
            // Flip hoặc Random angle
            (std::uniform_real_distribution<double>(0,1)(rng) < 0.6) ? (items[idx].angle + PI) : std::uniform_real_distribution<double>(0, 2*PI)(rng)
        );
    }
}

void SparrowSolver::explore() {
    auto start_time = std::chrono::steady_clock::now();
    auto deadline = start_time + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(config.TLx));
    
    if (!separate(config.Kx, config.Nx * 2, deadline)) {
        if (std::chrono::steady_clock::now() > deadline) return;
        config.container_size *= 1.1;
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
        for(auto& r : weights) std::fill(r.begin(), r.end(), 1.0); // Reset Master weights

        // Tăng Nx khi fail để kiên nhẫn hơn
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
            current_Rx *= 0.5; // Giảm bước nhảy
            scrambleItems(); 
            std::cout << "   -> Fail. Revert & Reduce Rx to " << current_Rx << std::endl;
        }
        iter++;
    }
    items = s_best;
    config.container_size = best_size;
}

void SparrowSolver::compress(double ratio) {
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
    std::cout << ">>> Solver Started. Items: " << items.size() << std::endl;
    int attempts = 0;
    while (attempts < 10) {
        auto deadline = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(5));
        if (separate(config.Kx, config.Nx * 2, deadline)) break;
        config.container_size *= 1.1; initializePlacement(); attempts++;
    }
    if (attempts >= 10) {
        std::cout << "Error: Could not find initial solution." << std::endl;
        return;
    }

    double bestSize = config.container_size;
    std::vector<CompositeShape> bestItems = items;
    
    for (int iter = 0; iter < config.max_restart; ++iter) {
        std::cout << "\n=== Major Iteration " << iter << " ===" << std::endl;
        descent();
        if (config.container_size < bestSize) { bestSize = config.container_size; bestItems = items; }
        items = bestItems; config.container_size = bestSize;
        
        explore();
        if (config.container_size < bestSize) { 
            bestSize = config.container_size; bestItems = items; 
            std::cout << "   [***] New Record: " << bestSize << std::endl;
        }
    }
    items = bestItems; config.container_size = bestSize;
    std::cout << ">>> Final Size: " << std::setprecision(10) << std::fixed << config.container_size << std::endl;
}