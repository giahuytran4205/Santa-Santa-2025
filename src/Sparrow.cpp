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
#include <numeric>

// ==========================================
// CONFIG LOADER (Giữ nguyên của bạn)
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
// SEPARATOR WORKER
// ==========================================

SeparatorWorker::SeparatorWorker(int _id, int seed, const std::vector<CompositeShape>& initItems, size_t n) 
    : id(_id), items(initItems), rng(seed) {
    weights.resize(n, std::vector<double>(n, 1.0));
    current_energy = 0.0;
    candidate_buffer.reserve(4096); 
}

void SeparatorWorker::load(const std::vector<CompositeShape>& masterItems, const std::vector<std::vector<double>>& masterWeights) {
    items = masterItems;
    weights = masterWeights;
}

double SeparatorWorker::randomDouble(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

// --- CỐT LÕI MỚI: Đánh giá bằng Circle Physics (Nhanh & Mượt) ---
double SeparatorWorker::evaluateSample(int itemIdx, Vec2 pos, double angle, const SparrowConfig& config) {
    // 1. Transform circles của item hiện tại ra World Space (On-the-fly)
    // Tận dụng struct SoA để tối ưu Cache
    const auto& localCircles = items[itemIdx].circles; 
    CirclesSoA transformedCircles;
    
    double c = std::cos(angle);
    double s = std::sin(angle);
    
    // Auto-vectorization loop
    for(size_t i = 0; i < localCircles.size(); ++i) {
        double lx = localCircles.x[i];
        double ly = localCircles.y[i];
        transformedCircles.add(
            lx * c - ly * s + pos.x,
            lx * s + ly * c + pos.y,
            localCircles.r[i]
        );
    }

    double totalCost = 0.0;
    double halfSize = config.container_size / 2.0;

    // 2. Wall Penalty (Dùng Circles để check biên cực nhanh)
    for(size_t i=0; i<transformedCircles.size(); ++i) {
        double x = transformedCircles.x[i];
        double y = transformedCircles.y[i];
        double r = transformedCircles.r[i];
        
        double pen = 0.0;
        if (x - r < -halfSize) pen += (-halfSize - (x - r));
        if (x + r > halfSize)  pen += ((x + r) - halfSize);
        if (y - r < -halfSize) pen += (-halfSize - (y - r));
        if (y + r > halfSize)  pen += ((y + r) - halfSize);
        
        if (pen > 0) {
            // Phạt lũy thừa để ép vật thể vào trong mạnh mẽ
            totalCost += (pen * 50.0) + (pen * pen * 500.0);
        }
    }

    // 3. Collision (Dùng evaluate_circles_fast thay cho SAT Polygon)
    for (size_t j = 0; j < items.size(); ++j) {
        if ((int)j == itemIdx) continue;
        
        // Transform item[j]
        CirclesSoA jWorld;
        double cj = std::cos(items[j].angle);
        double sj = std::sin(items[j].angle);
        const auto& jLocal = items[j].circles;
        
        for(size_t k=0; k<jLocal.size(); ++k) {
             jWorld.add(
                 jLocal.x[k] * cj - jLocal.y[k] * sj + items[j].pos.x,
                 jLocal.x[k] * sj + jLocal.y[k] * cj + items[j].pos.y,
                 jLocal.r[k]
             );
        }

        // Tính va chạm "mềm" (Soft Physics)
        double overlap = evaluate_circles_fast(transformedCircles, jWorld, 1e-15);
        if (overlap > 1e-20) {
            totalCost += overlap * weights[itemIdx][j];
        }
    }
    return totalCost;
}

// --- LOGIC TÌM KIẾM CỦA BẠN (Được giữ lại vì tốt, nhưng dùng hàm đánh giá mới) ---
void SeparatorWorker::searchPosition(int itemIdx, const SparrowConfig& config) {
    const CompositeShape& current = items[itemIdx];
    double current_e = evaluateSample(itemIdx, current.pos, current.angle, config);

    candidate_buffer.clear();
    
    double halfSize = config.container_size / 2.0;
    double foc_radius = config.container_size * 0.15; 

    // Sử dụng config.n_samples để chia pha tìm kiếm
    int n_div = (int)(config.n_samples * 0.4); 
    int n_foc = (int)(config.n_samples * 0.4); 
    int n_flip = config.n_samples - n_div - n_foc; 

    // Phase 1: Global Random
    for (int s = 0; s < n_div; ++s) {
        double x = randomDouble(-halfSize, halfSize);
        double y = randomDouble(-halfSize, halfSize);
        double rot = (randomDouble(0, 1) < 0.5) ? 0.0 : PI; 
        rot += randomDouble(-0.1, 0.1); 
        double e = evaluateSample(itemIdx, {x, y}, rot, config);
        candidate_buffer.emplace_back(e, Transform({x, y}, rot));
    }

    // Phase 2: Local Search
    for (int s = 0; s < n_foc; ++s) {
        double x = current.pos.x + randomDouble(-foc_radius, foc_radius);
        double y = current.pos.y + randomDouble(-foc_radius, foc_radius);
        double rot = current.angle + randomDouble(-0.1, 0.1);
        double e = evaluateSample(itemIdx, {x, y}, rot, config);
        candidate_buffer.emplace_back(e, Transform({x, y}, rot));
    }

    // Phase 3: Flip Strategy
    for (int s = 0; s < n_flip; ++s) {
        double x = current.pos.x + randomDouble(-foc_radius, foc_radius);
        double y = current.pos.y + randomDouble(-foc_radius, foc_radius);
        double rot = current.angle + PI + randomDouble(-0.1, 0.1);
        while (rot > 2*PI) rot -= 2*PI;
        double e = evaluateSample(itemIdx, {x, y}, rot, config);
        candidate_buffer.emplace_back(e, Transform({x, y}, rot));
    }

    std::sort(candidate_buffer.begin(), candidate_buffer.end(), 
        [](const auto& a, const auto& b){ return a.first < b.first; });

    // Local Refinement (Coordinate Descent) cho Top 3
    int K = 3; 
    Transform best_t(current.pos, current.angle);
    double best_e = current_e;

    for (int i = 0; i < std::min(K, (int)candidate_buffer.size()); ++i) {
        Transform t = candidate_buffer[i].second;
        double e = candidate_buffer[i].first;
        if (e > current_e * 1.5) continue; 

        double step = foc_radius * 0.5;
        double rot_step = 0.05; 
        
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
    // Sử dụng evaluateSample để nhất quán logic (chấp nhận tính dư 1 chút để code gọn)
    for (size_t i = 0; i < items.size(); ++i) {
        energy += evaluateSample((int)i, items[i].pos, items[i].angle, config);
    }
    return energy * 0.5; // Chia đôi vì overlap được tính 2 chiều
}

void SeparatorWorker::move_items(const SparrowConfig& config) {
    // 1. Chỉ di chuyển những item có năng lượng > 0 (đang va chạm hoặc lòi ra ngoài)
    std::vector<int> Ic;
    for (int i = 0; i < (int)items.size(); ++i) {
        // Threshold nhỏ để lọc
        if (evaluateSample(i, items[i].pos, items[i].angle, config) > 1e-6) {
            Ic.push_back(i);
        }
    }

    if (Ic.empty()) {
        current_energy = 0.0;
        return;
    }

    std::shuffle(Ic.begin(), Ic.end(), rng);

    // 2. Di chuyển
    for (int idx : Ic) {
        searchPosition(idx, config);
    }

    // 3. Tính lại energy
    current_energy = calculate_total_energy(config);
}


// ==========================================
// SPARROW SOLVER (MASTER)
// ==========================================

void SparrowSolver::saveState() {
    savedItems = items;
    savedContainerSize = config.container_size;
}

void SparrowSolver::loadState() {
    items = savedItems;
    config.container_size = savedContainerSize;
}

// --- TÍNH NĂNG MỚI: GENERATE SURROGATE CIRCLES ---
void SparrowSolver::generateSurrogateCircles() {
    std::cout << ">>> Generating Surrogate Circles..." << std::endl;
    for(auto& item : items) {
        item.generateSurrogateCircles(1.0); // 1.0 = Default quality
    }
}

// --- TÍNH NĂNG MỚI: LBF INITIALIZATION (Thay thế random grid) ---
// Dùng LBF để xếp hình chặt chẽ ngay từ đầu, giảm áp lực cho solver
void SparrowSolver::constructLBF() {
    std::cout << ">>> Running LBF Construction (Least Bad Fit)..." << std::endl;
    
    // 1. Sort: Ưu tiên vật lớn + dài
    std::vector<int> sortedIndices(items.size());
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(), [&](int a, int b) {
        return (items[a].convexHullArea * items[a].diameter) > 
               (items[b].convexHullArea * items[b].diameter);
    });

    // Reset ra xa
    for(auto& item : items) item.setTransform({999999, 999999}, 0);

    std::vector<int> placedIndices;
    placedIndices.reserve(items.size());
    double currentW = config.container_size;
    
    // Xếp lần lượt từng vật
    for (int itemIdx : sortedIndices) {
        CompositeShape& currentItem = items[itemIdx];
        double bestEnergy = std::numeric_limits<double>::max();
        Transform bestTransform = { {0,0}, 0 };
        int samples = 1000; 

        // Parallel sampling để tìm vị trí tốt nhất cho vật hiện tại
        #pragma omp parallel
        {
            double localMinE = std::numeric_limits<double>::max();
            Transform localBestT;
            int tid = omp_get_thread_num();
            std::mt19937 local_rng(rng() + tid + itemIdx);
            
            #pragma omp for
            for (int s = 0; s < samples; ++s) {
                double x = std::uniform_real_distribution<double>(-currentW/2, currentW/2)(local_rng);
                double y = std::uniform_real_distribution<double>(-currentW/2, currentW/2)(local_rng);
                double rot = (std::uniform_real_distribution<double>(0,1)(local_rng) < 0.7) ? 
                             0.0 : std::uniform_real_distribution<double>(0, 2*PI)(local_rng);

                // Quick Evaluation (Circle-based inline)
                CirclesSoA transCircles;
                double c = std::cos(rot), s_angle = std::sin(rot);
                const auto& orig = currentItem.circles;
                for(size_t k=0; k<orig.size(); ++k) {
                    transCircles.add(
                        orig.x[k]*c - orig.y[k]*s_angle + x,
                        orig.x[k]*s_angle + orig.y[k]*c + y,
                        orig.r[k]
                    );
                }
                
                double e = 0.0;
                // Wall check
                for(size_t k=0; k<transCircles.size(); ++k) {
                     double cx = transCircles.x[k], cy = transCircles.y[k], cr = transCircles.r[k];
                     if(cx-cr < -currentW/2 || cx+cr > currentW/2 || 
                        cy-cr < -currentW/2 || cy+cr > currentW/2) e += 10000.0;
                }
                
                // Collision check với các vật ĐÃ ĐẶT
                for (int pIdx : placedIndices) {
                    const auto& pItem = items[pIdx];
                    CirclesSoA pWorld; 
                    double pc = std::cos(pItem.angle), ps = std::sin(pItem.angle);
                    for(size_t k=0; k<pItem.circles.size(); ++k) {
                        pWorld.add(
                            pItem.circles.x[k]*pc - pItem.circles.y[k]*ps + pItem.pos.x,
                            pItem.circles.x[k]*ps + pItem.circles.y[k]*pc + pItem.pos.y,
                            pItem.circles.r[k]
                        );
                    }
                    e += evaluate_circles_fast(transCircles, pWorld, 1e-9);
                    if (e > localMinE) break; // Early exit optimization
                }

                if (e < localMinE) {
                    localMinE = e;
                    localBestT = { {x,y}, rot };
                }
            }
            #pragma omp critical
            {
                if (localMinE < bestEnergy) {
                    bestEnergy = localMinE;
                    bestTransform = localBestT;
                }
            }
        }
        
        currentItem.setTransform(bestTransform.pos, bestTransform.angle);
        placedIndices.push_back(itemIdx);
        
        // Nếu không nhét vừa -> Mở rộng container (Feature quan trọng của LBF)
        if (bestEnergy > 0.1) {
            currentW *= 1.1;
        }
    }
    // Cập nhật container size mới sau khi xếp xong
    config.container_size = currentW;
    std::cout << ">>> LBF Done. Initial Size: " << config.container_size << std::endl;
}

SparrowSolver::SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg) 
    : items(inputItems), config(cfg) {

    if (config.n_threads <= 0) config.n_threads = omp_get_num_procs();
    omp_set_num_threads(config.n_threads);

    std::random_device rd;
    rng = std::mt19937(rd());

    // 1. Generate Surrogate Circles (Bắt buộc cho engine mới)
    generateSurrogateCircles();

    // 2. Init Weights
    size_t n = items.size();
    weights.resize(n, std::vector<double>(n, 1.0));
    
    // 3. KHỞI TẠO BẰNG LBF (Thay vì Grid ngẫu nhiên)
    constructLBF();

    // 4. Init Persistent Workers
    workers.reserve(config.n_threads);
    for (int i = 0; i < config.n_threads; ++i) {
        workers.emplace_back(i, rd() + i, items, n);
    }
    
    std::cout << ">>> Solver Initialized with LBF & Circle Physics." << std::endl;
}

double SparrowSolver::getTotalEnergy() const {
    double energy = 0.0;
    // Tái sử dụng logic Circle Physics
    double halfSize = config.container_size / 2.0;

    for (size_t i = 0; i < items.size(); ++i) {
        // Transform on the fly
        CirclesSoA circlesI;
        double c = std::cos(items[i].angle), s = std::sin(items[i].angle);
        for(size_t k=0; k<items[i].circles.size(); ++k) {
            circlesI.add(items[i].circles.x[k]*c - items[i].circles.y[k]*s + items[i].pos.x,
                         items[i].circles.x[k]*s + items[i].circles.y[k]*c + items[i].pos.y,
                         items[i].circles.r[k]);
        }
        
        // Wall
        for(size_t k=0; k<circlesI.size(); ++k) {
             double x = circlesI.x[k], y = circlesI.y[k], r = circlesI.r[k];
             double pen = 0;
             if (x-r < -halfSize) pen += (-halfSize - (x-r));
             if (x+r > halfSize) pen += ((x+r) - halfSize);
             if (y-r < -halfSize) pen += (-halfSize - (y-r));
             if (y+r > halfSize) pen += ((y+r) - halfSize);
             if (pen>0) energy += pen*50 + pen*pen*500;
        }

        // Collision
        for (size_t j = i + 1; j < items.size(); ++j) {
            CirclesSoA circlesJ;
            double cj = std::cos(items[j].angle), sj = std::sin(items[j].angle);
            for(size_t k=0; k<items[j].circles.size(); ++k) {
                circlesJ.add(items[j].circles.x[k]*cj - items[j].circles.y[k]*sj + items[j].pos.x,
                             items[j].circles.x[k]*sj + items[j].circles.y[k]*cj + items[j].pos.y,
                             items[j].circles.r[k]);
            }
            energy += evaluate_circles_fast(circlesI, circlesJ, 1e-15);
        }
    }
    return energy;
}

void SparrowSolver::updateMasterWeights() {
    // Đơn giản hóa: reset trọng số về 1.0 hoặc giữ nguyên
    // Với circle physics mượt mà, GLS (Guided Local Search) ít quan trọng hơn so với SAT
    // Có thể implement sau nếu cần thiết.
}

// === ALGORITHM 9: SEPARATE ===
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

            // 1. SYNC & PARALLEL EXECUTION
            #pragma omp parallel for
            for (int i = 0; i < (int)workers.size(); ++i) {
                workers[i].load(items, weights);
                workers[i].move_items(config); // Dùng config.n_samples, etc.
            }

            // 2. GATHER RESULT
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
                if (min_worker_energy < best_global_energy) {
                    items = workers[best_worker_idx].items;
                    best_global_energy = min_worker_energy;
                    improved_in_try = true;
                }
            }

            if (best_global_energy <= 1e-5) return true;

            // Stagnation logic
            if (min_worker_energy < best_global_energy * 0.999) { 
                n = 0; 
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
    int count = std::max(1, (int)(items.size() * 0.3)); 
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
            (std::uniform_real_distribution<double>(0,1)(rng) < 0.6) ? (items[idx].angle + PI) : std::uniform_real_distribution<double>(0, 2*PI)(rng)
        );
    }
}

// Logic Explore/Compress/Descent được giữ nguyên cấu trúc
// nhưng giờ nó gọi separate() với engine circle physics mạnh mẽ bên dưới
void SparrowSolver::explore() {
    auto start_time = std::chrono::steady_clock::now();
    auto deadline = start_time + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(config.TLx));
    
    if (!separate(config.Kx, config.Nx * 2, deadline)) {
        if (std::chrono::steady_clock::now() > deadline) return;
        config.container_size *= 1.1;
        constructLBF(); // QUAN TRỌNG: Dùng LBF để reset khi tắc
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
        config.container_size *= 1.1; 
        constructLBF(); // QUAN TRỌNG: Sử dụng LBF để xếp lại nếu thất bại
        attempts++;
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