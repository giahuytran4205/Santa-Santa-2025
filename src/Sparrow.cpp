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

// --- Helper Functions ---

// Thread-Safe Random Helper
double SparrowSolver::randomDouble(double min, double max) {
    static thread_local std::mt19937* local_rng = nullptr;
    if (!local_rng) {
        // Seed based on time and thread ID to ensure independence
        unsigned int seed = (unsigned int)(std::chrono::system_clock::now().time_since_epoch().count() + omp_get_thread_num());
        local_rng = new std::mt19937(seed);
    }
    std::uniform_real_distribution<double> dist(min, max);
    return dist(*local_rng);
}

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
            
            // Chessboard heuristic
            double rot = ((r + c) % 2 == 0) ? 0.0 : PI;
            
            rot += randomDouble(-0.1, 0.1); 
            double jitter = spacing * 0.1;
            x += randomDouble(-jitter, jitter);
            y += randomDouble(-jitter, jitter);

            items[itemIdx].setTransform({x, y}, rot);
            idx++;
        }
    }
    std::cout << "Initialized placement (Chessboard heuristic), initial energy: " << getTotalEnergy() << std::endl;
}

SparrowSolver::SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg) 
    : items(inputItems), config(cfg) {
        
    if (config.n_threads > 0) {
        omp_set_num_threads(config.n_threads);
    } else {
        omp_set_num_threads(omp_get_num_procs());
    }

    #pragma omp parallel
    {
        #pragma omp single
        {
            std::cout << ">>> OpenMP Optimized Solver Running with " << omp_get_num_threads() << " threads." << std::endl;
        }
    }
    
    std::random_device rd;
    rng = std::mt19937(rd());

    size_t n = items.size();
    weights.resize(n, std::vector<double>(n, 1.0));

    // Precompute poles
    for (auto& item : items) {
        item.precomputeAllPoles(0.01); 
    }

    initializePlacement();
}

// Algorithm 7: Evaluate Sample
double SparrowSolver::evaluateSample(int itemIdx, Vec2 pos, double angle, const std::vector<CompositeShape>& local_items, const std::vector<std::vector<double>>& local_weights) const {
    CompositeShape tempItem = local_items[itemIdx];
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
        totalCost += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
    }

    // 2. Item Overlap
    for (size_t j = 0; j < local_items.size(); ++j) {
        if (static_cast<int>(j) == itemIdx) continue;
        
        if (tempItem.totalAABB.overlaps(local_items[j].totalAABB)) {
             double overlap = quantify_collision(tempItem, local_items[j]);
             if (overlap > 1e-9) {
                 totalCost += overlap * local_weights[itemIdx][j];
             }
        }
    }
    return totalCost;
}

// --- OPTIMIZED SEARCH: Sequential logic (called by parallel workers) ---
void SparrowSolver::searchPosition(int itemIdx, std::vector<CompositeShape>& local_items, std::vector<std::vector<double>>& local_weights) {
    const CompositeShape& current = local_items[itemIdx];
    double current_e = evaluateSample(itemIdx, current.pos, current.angle, local_items, local_weights);

    static thread_local std::vector<std::pair<double, Transform>> candidates;
    candidates.clear(); // Xóa dữ liệu cũ, nhưng giữ nguyên capacity bộ nhớ
    candidates.reserve(config.n_samples);

    double halfSize = config.container_size / 2.0;
    double foc_radius = config.container_size * 0.15;
    
    int n_div = (int)(config.n_samples * 0.4);
    int n_foc = (int)(config.n_samples * 0.4);
    int n_flip = config.n_samples - n_div - n_foc;

    // Phase 1: Global
    for (int s = 0; s < n_div; ++s) {
        double x = randomDouble(-halfSize, halfSize);
        double y = randomDouble(-halfSize, halfSize);
        double rot = ((int)(randomDouble(0, 100)) % 2 == 0) ? 0.0 : PI;
        rot += randomDouble(-0.01, 0.01);
        double e = evaluateSample(itemIdx, {x, y}, rot, local_items, local_weights);
        candidates.emplace_back(e, Transform({x, y}, rot));
    }

    // Phase 2: Local
    for (int s = 0; s < n_foc; ++s) {
        double dx = randomDouble(-foc_radius, foc_radius);
        double dy = randomDouble(-foc_radius, foc_radius);
        double rot = current.angle + randomDouble(-0.02, 0.02);
        double e = evaluateSample(itemIdx, current.pos + Vec2{dx, dy}, rot, local_items, local_weights);
        candidates.emplace_back(e, Transform(current.pos + Vec2{dx, dy}, rot));
    }

    // Phase 3: Flip
    for (int s = 0; s < n_flip; ++s) {
        double dx = randomDouble(-foc_radius, foc_radius);
        double dy = randomDouble(-foc_radius, foc_radius);
        double rot = current.angle + PI + randomDouble(-0.01, 0.01);
        while (rot > 2 * PI) rot -= 2 * PI;
        double e = evaluateSample(itemIdx, current.pos + Vec2{dx, dy}, rot, local_items, local_weights);
        candidates.emplace_back(e, Transform(current.pos + Vec2{dx, dy}, rot));
    }

    // Sort best candidates
    std::sort(candidates.begin(), candidates.end());

    // Coordinate Descent on top K candidates
    int K = 3; 
    Transform best_t(current.pos, current.angle);
    double best_e = current_e;

    for (int i = 0; i < std::min(K, (int)candidates.size()); ++i) {
        Transform t = candidates[i].second;
        double e = candidates[i].first;
        if (e > current_e * 1.5) continue;

        double step = foc_radius * 0.5;
        double rot_step = 0.05; 
        
        for (int iter = 0; iter < 20; ++iter) { 
             bool improved = false;
             // X
             for (double d : {-step, step}) {
                 double ne = evaluateSample(itemIdx, t.pos + Vec2{d,0}, t.angle, local_items, local_weights);
                 if (ne < e) { t.pos.x += d; e = ne; improved = true; }
             }
             // Y
             for (double d : {-step, step}) {
                 double ne = evaluateSample(itemIdx, t.pos + Vec2{0,d}, t.angle, local_items, local_weights);
                 if (ne < e) { t.pos.y += d; e = ne; improved = true; }
             }
             // Angle
             for (double d : {-rot_step, rot_step}) {
                 double ne = evaluateSample(itemIdx, t.pos, t.angle + d, local_items, local_weights);
                 if (ne < e) { t.angle += d; e = ne; improved = true; }
             }
             if (!improved) { step *= 0.5; rot_step *= 0.5; }
             if (step < 1e-4) break;
        }

        if (e < best_e) { best_e = e; best_t = t; }
    }

    if (best_e < current_e) {
        local_items[itemIdx].setTransform(best_t.pos, best_t.angle);
    }
}

double SparrowSolver::getTotalEnergy() {
    return calculate_energy(items);
}

double SparrowSolver::calculate_energy(const std::vector<CompositeShape>& local_items) const {
    double energy = 0.0;
    
    // Item-Item Overlap
    for (size_t a = 0; a < local_items.size(); ++a) {
        for (size_t b = a + 1; b < local_items.size(); ++b) {
            if (local_items[a].totalAABB.overlaps(local_items[b].totalAABB)) {
                energy += quantify_collision(local_items[a], local_items[b]);
            }
        }
    }
    
    // Wall Penalty
    double halfSize = config.container_size / 2.0;
    for (const auto& item : local_items) {
        const AABB& aabb = item.totalAABB;
        double outOfBounds = 0.0;
        if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
        if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
        if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
        if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
        
        if (outOfBounds > 0) {
            energy += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
        }
    }
    return energy;
}

void SparrowSolver::updateWeights(std::vector<std::vector<double>>& local_weights) {
    // This is a placeholder for the class method, but we use the helper below in parallel.
}

// Helper to update weights inside worker
void updateLocalWeights(std::vector<std::vector<double>>& w, const std::vector<CompositeShape>& it, const SparrowConfig& cfg) {
    double e_max = 0.0;
    size_t n = it.size();
    for (size_t a = 0; a < n; ++a) {
        for (size_t b = a + 1; b < n; ++b) {
            if (it[a].totalAABB.overlaps(it[b].totalAABB)) {
                double e = quantify_collision(it[a], it[b]);
                if (e > e_max) e_max = e;
            }
        }
    }
    
    if (e_max < 1e-9) return;

    for (size_t a = 0; a < n; ++a) {
        for (size_t b = a + 1; b < n; ++b) {
            double e = 0.0;
            if (it[a].totalAABB.overlaps(it[b].totalAABB)) {
                e = quantify_collision(it[a], it[b]);
            }
            
            double m;
            if (e > 1e-9) m = cfg.Ml + (cfg.Mu - cfg.Ml) * (e / e_max);
            else m = cfg.Md;
            
            w[a][b] = std::max(1.0, w[a][b] * m);
            w[b][a] = w[a][b]; 
        }
    }
}

void SparrowSolver::perform_move_items(std::vector<CompositeShape>& local_items, std::vector<std::vector<double>>& local_weights) {
    std::vector<int> Ic;
    double halfSize = config.container_size / 2.0;
    
    // Identify Colliding Items
    for (int i = 0; i < local_items.size(); ++i) {
        bool colliding = false;
        const AABB& aabb = local_items[i].totalAABB;
        
        if (aabb.min.x < -halfSize || aabb.max.x > halfSize || aabb.min.y < -halfSize || aabb.max.y > halfSize) {
            colliding = true;
        } else {
            for (int j = 0; j < local_items.size(); ++j) {
                if (i == j) continue;
                if (local_items[i].totalAABB.overlaps(local_items[j].totalAABB)) {
                    if (quantify_collision(local_items[i], local_items[j]) > 1e-6) {
                        colliding = true;
                        break;
                    }
                }
            }
        }
        if (colliding) Ic.push_back(i);
    }

    if (Ic.empty()) return;

    // --- FIX: Use correct thread-local RNG ---
    static thread_local std::mt19937* local_rng_ptr = nullptr;
    if(!local_rng_ptr) {
         unsigned int seed = (unsigned int)(std::chrono::system_clock::now().time_since_epoch().count() + omp_get_thread_num());
         local_rng_ptr = new std::mt19937(seed);
    }
    std::shuffle(Ic.begin(), Ic.end(), *local_rng_ptr);

    for (int i : Ic) {
        searchPosition(i, local_items, local_weights);
    }
}

// --- CORE PARALLEL OPTIMIZATION ---
bool SparrowSolver::separate(int kmax, int nmax, double time_limit_sec) {
    auto start = std::chrono::steady_clock::now();
    
    // Initial best state
    double best_global_energy = getTotalEnergy();
    if (best_global_energy <= 1e-4) return true;

    // Flag to stop early if feasible solution found
    bool solution_found = false;

    // Start Parallel Region
    // Each thread becomes an independent worker (Walker)
    #pragma omp parallel
    {
        // 1. Thread-Local State Initialization (Copy from Master)
        static thread_local std::vector<CompositeShape> local_items;
        static thread_local std::vector<std::vector<double>> local_weights;

        // Phép gán này sẽ tái sử dụng capacity cũ, không cần malloc lại
        local_items = items;
        local_weights = weights;
        
        // 2. Worker Loop
        int k = 0;
        while (k < kmax && !solution_found) {
            
            int n = 0;
            bool improved_in_try = false;
            double local_best_energy = calculate_energy(local_items);

            while (n < nmax) {
                // Check time limit periodically
                if (n % 10 == 0) {
                    if (solution_found) break; // Break inner loop
                    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count();
                    if (elapsed > time_limit_sec) break;
                }

                // Move items (Sequential logic within worker)
                perform_move_items(local_items, local_weights);
                
                // Update GLS weights
                updateLocalWeights(local_weights, local_items, config);

                double current_e = calculate_energy(local_items);

                // Update Local Best
                if (current_e < local_best_energy) {
                    local_best_energy = current_e;
                    improved_in_try = true;
                }

                // Check against Global Best
                // Use double-checked locking to avoid mutex overhead
                if (current_e < best_global_energy) {
                    #pragma omp critical
                    {
                        if (current_e < best_global_energy) {
                            best_global_energy = current_e;
                            items = local_items; // Update Master State
                            weights = local_weights;
                            std::cout << "   [Worker " << omp_get_thread_num() << "] New Best Energy: " << best_global_energy << std::endl;
                        }
                        if (best_global_energy <= 1e-4) {
                            solution_found = true;
                        }
                    }
                }

                if (solution_found) break;

                // Stagnation counter logic
                if (current_e >= local_best_energy) n++;
                else n = 0;
            }
            
            if (solution_found) break;

            k++;
            if (improved_in_try) k = 0;
        }
    } // End Parallel

    return (best_global_energy <= 1e-4);
}

void SparrowSolver::move_items_multi() {
    // Deprecated / No-op as logic is moved to separate()
}

void SparrowSolver::scrambleItems() {
    int count = std::max(1, (int)(items.size() * 0.2));
    std::vector<int> indices(items.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    double halfSize = config.container_size / 2.0;

    for (int i = 0; i < count; ++i) {
        int idx = indices[i];
        Vec2 newPos;
        newPos.x = randomDouble(-halfSize, halfSize);
        newPos.y = randomDouble(-halfSize, halfSize);
        
        double currentAngle = items[idx].angle;
        double newAngle;
        if (randomDouble(0, 1) < 0.5) newAngle = currentAngle + PI;
        else newAngle = randomDouble(0, 2.0 * PI);
        newAngle += randomDouble(-0.1, 0.1);

        items[idx].setTransform(newPos, newAngle);
    }
    std::cout << "   [Scramble] Perturbed " << count << " items." << std::endl;
}

void SparrowSolver::explore() {
    auto start = std::chrono::steady_clock::now();
    
    if (!separate(config.Kx, config.Nx * 2, config.TLx)) {
        config.container_size *= 1.05;
        initializePlacement(); 
        separate(config.Kx, config.Nx * 2, config.TLx);
    }

    double best_size = config.container_size;
    std::vector<CompositeShape> s_best = items;

    int iter = 0;
    int fail_streak = 0;
    double current_Rx = config.Rx; 

    while (iter < config.max_outer_loops) {
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count();
        if (elapsed > config.TLx) break;
        
        if (current_Rx < 1e-5) {
             current_Rx = config.Rx; 
             scrambleItems(); 
        }

        std::vector<CompositeShape> s_prev_feasible = items;
        double size_prev = config.container_size;

        config.container_size *= (1.0 - current_Rx);

        std::cout << "Explore [" << iter << "] | Rx: " << current_Rx 
                  << " | Target: " << config.container_size 
                  << " | Streak: " << fail_streak << std::endl;

        double scale = config.container_size / size_prev;
        for(auto& it : items) it.setTransform(it.pos * scale, it.angle);
        
        for(auto& r : weights) std::fill(r.begin(), r.end(), 1.0);

        int effective_Nx = config.Nx + (fail_streak * 50);
        bool feasible = separate(config.Kx, effective_Nx, config.TLx - elapsed);

        if (feasible) {
            best_size = config.container_size;
            s_best = items;
            fail_streak = 0;
        } else {
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
    auto start = std::chrono::steady_clock::now();
    std::vector<CompositeShape> s_star = items;
    double best_size = config.container_size;
    
    int loop_iter = 0;
    while (loop_iter < config.max_outer_loops) {
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count();
        if (elapsed > config.TLc) break;

        double progress = (double)elapsed / config.TLc;
        double r = config.Rs_c + (config.Re_c - config.Rs_c) * progress;

        items = s_star;
        double old_size = config.container_size;
        config.container_size *= (1.0 - r);

        for(auto& row : weights) std::fill(row.begin(), row.end(), 1.0);

        bool feasible = separate(config.Kc, config.Nc, config.TLc - elapsed);

        if (feasible) {
            s_star = items;
            best_size = config.container_size;
             std::cout << "   Compress success -> " << best_size << std::endl;
        } else {
            config.container_size = old_size;
             std::cout << "   Compress fail, revert." << std::endl;
        }
        loop_iter++;
    }
    items = s_star;
    config.container_size = best_size;
}

void SparrowSolver::descent() {
    compress(0.0);
}

void SparrowSolver::solve() {
    std::cout << ">>> Solver Started. Items: " << items.size() << std::endl;
    
    int attempts = 0;
    while (!separate(config.Kx, config.Nx * 2, config.TLx) && attempts < 10) {
        std::cout << "[Warn] Initial config infeasible. Expanding..." << std::endl;
        config.container_size *= 1.1;
        initializePlacement();
        attempts++;
    }
    
    if (attempts >= 10) {
        std::cout << "[Error] Could not find initial feasible solution." << std::endl;
        return;
    }

    std::cout << ">>> Initial Feasible Found: " << config.container_size << std::endl;

    double bestSize = config.container_size;
    std::vector<CompositeShape> bestItems = items;
    
    int maxRestarts = config.max_restart;
    for (int iter = 0; iter < maxRestarts; ++iter) {
        std::cout << "\n=== Major Iteration " << iter << " ===" << std::endl;
        descent();
        
        if (config.container_size < bestSize) {
            bestSize = config.container_size;
            bestItems = items;
            std::cout << "   [***] New Record: " << bestSize << std::endl;
        }

        items = bestItems;
        config.container_size = bestSize;
        explore();
        
        if (config.container_size < bestSize) {
            bestSize = config.container_size;
            bestItems = items;
            std::cout << "   [***] Explore found better: " << bestSize << std::endl;
        }
    }

    items = bestItems;
    config.container_size = bestSize;
    std::cout << ">>> Final Size: " << std::setprecision(10) << std::fixed << config.container_size << std::endl;
}