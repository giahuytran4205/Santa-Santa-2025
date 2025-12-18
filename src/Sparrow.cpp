// Sparrow.cpp
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
        std::random_device rd;
        local_rng = new std::mt19937(rd() + omp_get_thread_num());
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

// Khởi tạo thông minh hơn: Xếp lưới nhưng xen kẽ hướng (Up/Down)
void SparrowSolver::initializePlacement() {
    double halfSize = config.container_size / 2.0;
    int n = items.size();
    int rows = std::ceil(std::sqrt(n));
    double spacing = config.container_size / rows;
    int idx = 0;
    
    // Shuffle thứ tự để tránh bias
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < rows && idx < n; ++c) {
            int itemIdx = indices[idx];
            double x = -halfSize + spacing * (c + 0.5);
            double y = -halfSize + spacing * (r + 0.5);
            
            // Heuristic: Lưới bàn cờ (Chessboard pattern) cho hướng xoay
            // Một cây hướng lên (0), cây bên cạnh hướng xuống (PI)
            double rot = ((r + c) % 2 == 0) ? 0.0 : PI;
            
            // Thêm nhiễu nhẹ để không quá cứng nhắc
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

// Constructor
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
            std::cout << ">>> OpenMP Running with " << omp_get_num_threads() << " threads." << std::endl;
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

    // 1. Wall Penalty (Soft Constraint)
    // Giảm penalty xuống để gradient mượt hơn (10000 -> 100 * overlap)
    const AABB& aabb = tempItem.totalAABB;
    double outOfBounds = 0.0;
    if (aabb.min.x < -halfSize) outOfBounds += (-halfSize - aabb.min.x);
    if (aabb.max.x > halfSize)  outOfBounds += (aabb.max.x - halfSize);
    if (aabb.min.y < -halfSize) outOfBounds += (-halfSize - aabb.min.y);
    if (aabb.max.y > halfSize)  outOfBounds += (aabb.max.y - halfSize);
    
    // Penalty lũy thừa để càng ra xa càng bị phạt nặng, nhưng ở gần biên thì nhẹ nhàng
    if (outOfBounds > 0) {
        totalCost += (outOfBounds * 100.0) + (outOfBounds * outOfBounds * 1000.0);
    }

    // 2. Item Overlap (Weighted by GLS)
    for (size_t j = 0; j < local_items.size(); ++j) {
        if (static_cast<int>(j) == itemIdx) continue;
        
        // Chỉ check broadphase AABB trước cho nhanh
        if (tempItem.totalAABB.overlaps(local_items[j].totalAABB)) {
             double overlap = quantify_collision(tempItem, local_items[j]);
             if (overlap > 1e-9) {
                 totalCost += overlap * local_weights[itemIdx][j];
             }
        }
    }
    return totalCost;
}

// Algorithm 6: Search Position (Đã cải tiến cho Homogeneous/Symmetry)
void SparrowSolver::searchPosition(int itemIdx, std::vector<CompositeShape>& local_items, std::vector<std::vector<double>>& local_weights) {
    const CompositeShape& current = local_items[itemIdx];
    Transform current_t(current.pos, current.angle);
    double current_e = evaluateSample(itemIdx, current_t.pos, current_t.angle, local_items, local_weights);

    std::vector<std::pair<double, Transform>> samples;
    samples.reserve(config.n_samples);

    double halfSize = config.container_size / 2.0;
    double foc_radius = config.container_size * 0.15; // Tăng vùng tìm kiếm cục bộ lên chút

    // 1. Global Random Samples (T_div) - 40%
    int n_div = (int)(config.n_samples * 0.4);
    for (int s = 0; s < n_div; ++s) {
        double x = randomDouble(-halfSize, halfSize);
        double y = randomDouble(-halfSize, halfSize);
        
        // Bias Rotation: Ưu tiên 0 và PI (Up/Down) cho cây
        double r = randomDouble(0, 1);
        double rot;
        if (r < 0.7) rot = (randomDouble(0, 1) < 0.5) ? 0.0 : PI; 
        else rot = randomDouble(0, 2.0 * PI);
        
        // Thêm nhiễu góc
        rot += randomDouble(-0.15, 0.15); 

        double e = evaluateSample(itemIdx, {x, y}, rot, local_items, local_weights);
        samples.emplace_back(e, Transform({x, y}, rot));
    }

    // 2. Local Samples (T_foc) - 40%
    int n_foc = (int)(config.n_samples * 0.4);
    for (int s = 0; s < n_foc; ++s) {
        double dx = randomDouble(-foc_radius, foc_radius);
        double dy = randomDouble(-foc_radius, foc_radius);
        double x = current.pos.x + dx;
        double y = current.pos.y + dy;
        
        // Giữ góc hiện tại + nhiễu
        double drot = randomDouble(-0.2, 0.2); 
        double rot = current.angle + drot;

        double e = evaluateSample(itemIdx, {x, y}, rot, local_items, local_weights);
        samples.emplace_back(e, Transform({x, y}, rot));
    }

    // 3. FLIP SAMPLES (Quan trọng nhất cho đối xứng trục) - 20%
    // Thử lật ngược vật thể tại vị trí hiện tại (hoặc lân cận)
    int n_flip = config.n_samples - n_div - n_foc;
    for (int s = 0; s < n_flip; ++s) {
        double dx = randomDouble(-foc_radius*0.5, foc_radius*0.5);
        double dy = randomDouble(-foc_radius*0.5, foc_radius*0.5);
        double x = current.pos.x + dx;
        double y = current.pos.y + dy;

        // LẬT 180 độ so với góc hiện tại
        double rot = current.angle + PI + randomDouble(-0.1, 0.1); 

        double e = evaluateSample(itemIdx, {x, y}, rot, local_items, local_weights);
        samples.emplace_back(e, Transform({x, y}, rot));
    }

    // Sort chọn cái tốt nhất
    std::sort(samples.begin(), samples.end());

    // Refinement (Coordinate Descent) cho top K
    int K = 3; 
    Transform best_t = current_t;
    double best_e = current_e;

    for (int i = 0; i < std::min(K, (int)samples.size()); ++i) {
        Transform t = samples[i].second;
        double e = samples[i].first;

        // Nếu mẫu thử tệ hơn nhiều so với hiện tại thì bỏ qua (trừ khi đang early stage)
        if (e > current_e * 1.5) continue;

        // Adaptive Descent
        double step = foc_radius * 0.5;
        double rot_step = 0.1; // ~5.7 độ
        
        for (int iter = 0; iter < 10; ++iter) { 
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
            
            if (!improved) {
                step *= 0.5;
                rot_step *= 0.5;
            }
            if (step < 1e-4) break;
        }

        if (e < best_e) {
            best_e = e;
            best_t = t;
        }
    }

    // Apply best
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
            // Broadphase check
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
    double e_max = 0.0;
    // Tìm max overlap hiện tại
    for (size_t a = 0; a < items.size(); ++a) {
        for (size_t b = a + 1; b < items.size(); ++b) {
            if (items[a].totalAABB.overlaps(items[b].totalAABB)) {
                double e = quantify_collision(items[a], items[b]);
                if (e > e_max) e_max = e;
            }
        }
    }
    
    if (e_max < 1e-9) return;

    for (size_t a = 0; a < items.size(); ++a) {
        for (size_t b = a + 1; b < items.size(); ++b) {
            double e = 0.0;
            if (items[a].totalAABB.overlaps(items[b].totalAABB)) {
                e = quantify_collision(items[a], items[b]);
            }
            
            double m;
            if (e > 1e-9) {
                // Tăng trọng số cho cặp va chạm mạnh
                m = config.Ml + (config.Mu - config.Ml) * (e / e_max);
            } else {
                // Giảm dần trọng số cho cặp đã an toàn
                m = config.Md;
            }
            local_weights[a][b] = std::max(1.0, local_weights[a][b] * m);
            local_weights[b][a] = local_weights[a][b]; 
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
        
        // Wall check
        if (aabb.min.x < -halfSize || aabb.max.x > halfSize || aabb.min.y < -halfSize || aabb.max.y > halfSize) {
            colliding = true;
        }
        
        // Item check (nếu chưa chạm tường)
        if (!colliding) {
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

    std::random_device rd;
    std::mt19937 local_rng(rd());
    std::shuffle(Ic.begin(), Ic.end(), local_rng);

    for (int i : Ic) {
        searchPosition(i, local_items, local_weights);
    }
}

void SparrowSolver::move_items_multi() {
    // Biến chia sẻ (Shared) để lưu kết quả tốt nhất
    double best_e = std::numeric_limits<double>::max();
    
    // Lưu ý: Không khởi tạo best_items ở đây để tránh copy không cần thiết nếu không tìm thấy cái mới
    // Ta dùng cờ hoặc mutex để update sau.
    // Tuy nhiên, logic cũ dùng critical section là an toàn.
    
    // Biến tạm để lưu kết quả tốt nhất tìm được trong lần chạy này (để hạn chế critical section)
    std::vector<CompositeShape> current_best_items;
    std::vector<std::vector<double>> current_best_weights;
    bool found_improvement = false;

    #pragma omp parallel 
    {
        // --- TỐI ƯU HÓA BỘ NHỚ ---
        // Sử dụng static thread_local để biến này "sống" mãi theo từng luồng.
        // Nó chỉ được cấp phát 1 lần duy nhất khi luồng khởi tạo, và được dùng lại mãi mãi.
        static thread_local std::vector<CompositeShape> local_items;
        static thread_local std::vector<std::vector<double>> local_weights;

        // Gán dữ liệu mới vào vùng nhớ cũ (Reuse capacity)
        // std::vector::operator= sẽ tái sử dụng vùng nhớ đã cấp phát nếu đủ chỗ.
        // Điều này biến chi phí từ "malloc + copy" thành chỉ còn "copy" (nhanh hơn rất nhiều).
        local_items = items; 
        local_weights = weights;

        // --- LOGIC XỬ LÝ ---
        perform_move_items(local_items, local_weights);
        updateWeights(local_weights);

        double local_e = calculate_energy(local_items);
        
        // --- GOM KẾT QUẢ ---
        // Kiểm tra nhanh trước khi vào critical section để tránh nghẽn cổ chai
        if (local_e < best_e) {
            #pragma omp critical
            {
                // Kiểm tra lại lần nữa trong vùng an toàn (Double-checked locking pattern)
                if (local_e < best_e) {
                    best_e = local_e;
                    // Copy ra ngoài vùng shared
                    // Lưu ý: Việc copy này ít xảy ra hơn nhiều so với việc tính toán
                    current_best_items = local_items;
                    current_best_weights = local_weights;
                    found_improvement = true;
                }
            }
        }
    }

    // Chỉ cập nhật vào items chính nếu tìm thấy cải thiện
    if (found_improvement) {
        items = current_best_items;
        weights = current_best_weights;
    }
}

bool SparrowSolver::separate(int kmax, int nmax, double time_limit_sec) {
    auto start = std::chrono::steady_clock::now();
    double e_star = getTotalEnergy();
    std::vector<CompositeShape> s_star = items; 
    
    int k = 0;
    
    // Nếu đã feasible ngay từ đầu
    if (e_star <= 1e-4) return true;

    while (k < kmax) {
        // Reset về state tốt nhất mỗi lần thử lại (Strike system)
        items = s_star; 
        
        int n = 0;
        bool improved_in_try = false;

        while (n < nmax) {
            move_items_multi(); // Cập nhật items và weights
            double e = getTotalEnergy();
            
            if (e < e_star) {
                e_star = e;
                s_star = items;
                improved_in_try = true;
                // Nếu tìm thấy feasible thì return ngay
                if (e_star <= 1e-4) {
                    items = s_star;
                    return true;
                }
            }

            // Tăng counter nếu không cải thiện
            if (e >= e_star) {
                n++;
            } else {
                // Nếu có cải thiện, reset counter để đào sâu tiếp
                n = 0; 
            }

            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count();
            if (elapsed > time_limit_sec) {
                items = s_star;
                return (e_star <= 1e-4);
            }
        }
        
        // Hết nmax bước mà không feasible -> Strike
        k++;
        // Nếu trong lần thử vừa rồi có tiến bộ, reset strike để cho thêm cơ hội
        if (improved_in_try) k = 0; 
        
        std::cout << "   Separate step | Energy: " << e_star << " | Strike: " << k << "/" << kmax << std::endl;
    }

    items = s_star;
    return (e_star <= 1e-4);
}

// FIX: Hàm Scramble không Swap nữa mà Perturb/Flip
void SparrowSolver::scrambleItems() {
    // Chọn ngẫu nhiên 20% số lượng vật thể để gây nhiễu
    int count = std::max(1, (int)(items.size() * 0.2));
    std::vector<int> indices(items.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    double halfSize = config.container_size / 2.0;

    for (int i = 0; i < count; ++i) {
        int idx = indices[i];
        
        // 1. Dịch chuyển ngẫu nhiên
        Vec2 newPos;
        newPos.x = randomDouble(-halfSize, halfSize);
        newPos.y = randomDouble(-halfSize, halfSize);

        // 2. XOAY/LẬT (Quan trọng)
        // Thay đổi hoàn toàn trạng thái orientation để phá kẹt
        double currentAngle = items[idx].angle;
        double newAngle;
        
        if (randomDouble(0, 1) < 0.5) {
            newAngle = currentAngle + PI; // Flip 180
        } else {
            newAngle = randomDouble(0, 2.0 * PI); // Random mới
        }
        // Thêm nhiễu
        newAngle += randomDouble(-0.1, 0.1);

        items[idx].setTransform(newPos, newAngle);
    }
    std::cout << "   [Scramble] Perturbed " << count << " items." << std::endl;
}

void SparrowSolver::explore() {
    auto start = std::chrono::steady_clock::now();
    
    // Khởi đầu: Tìm feasible đầu tiên
    if (!separate(config.Kx, config.Nx * 2, config.TLx)) {
        config.container_size *= 1.05;
        initializePlacement(); 
        separate(config.Kx, config.Nx * 2, config.TLx);
    }

    double best_size = config.container_size;
    std::vector<CompositeShape> s_best = items;

    int iter = 0;
    int fail_streak = 0; // Đếm số lần thất bại liên tiếp
    
    // Sử dụng current_Rx riêng để có thể thay đổi linh hoạt
    double current_Rx = config.Rx; 

    while (iter < config.max_outer_loops) {
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count();
        if (elapsed > config.TLx) break;
        
        // Nếu Rx quá nhỏ (dưới 0.001%), ta coi như đã hội tụ, reset lại để explore vùng khác
        if (current_Rx < 1e-5) {
             current_Rx = config.Rx; // Reset step size
             // Có thể scramble nhẹ để đổi gió
             scrambleItems(); 
        }

        // Lưu trạng thái Feasible tốt nhất hiện tại
        std::vector<CompositeShape> s_prev_feasible = items;
        double size_prev = config.container_size;

        // --- ADAPTIVE SHRINK LOGIC ---
        // Thu nhỏ thùng chứa
        config.container_size *= (1.0 - current_Rx);

        std::cout << "Explore [" << iter << "] | Rx: " << current_Rx 
                  << " | Target: " << config.container_size 
                  << " | Streak: " << fail_streak << std::endl;

        // Scale vị trí các vật vào trong
        double scale = config.container_size / size_prev;
        for(auto& it : items) it.setTransform(it.pos * scale, it.angle);
        
        // Reset GLS weights
        for(auto& r : weights) std::fill(r.begin(), r.end(), 1.0);

        // Cố gắng giải
        // Nếu fail streak cao, tăng nỗ lực giải (Nx) lên để cố cứu vãn
        int effective_Nx = config.Nx + (fail_streak * 50);
        bool feasible = separate(config.Kx, effective_Nx, config.TLx - elapsed);

        if (feasible) {
            // THÀNH CÔNG
            best_size = config.container_size;
            s_best = items;
            
            // Nếu thành công, ta có thể giữ nguyên Rx hoặc tăng nhẹ để đi nhanh hơn (Momentum)
            fail_streak = 0;
            // current_Rx = std::min(config.Rx, current_Rx * 1.1); 
        } else {
            // THẤT BẠI (Infeasible)
            fail_streak++;
            
            // Revert về trạng thái tốt nhất trước đó
            items = s_prev_feasible;
            config.container_size = size_prev;
            
            // QUAN TRỌNG: Giảm Rx đi một nửa để lần sau thử bước nhỏ hơn
            current_Rx *= 0.5;
            
            // Đồng thời Scramble để thay đổi cấu hình, hy vọng cấu hình mới sẽ chịu được mức nén này
            // Nhưng chỉ Scramble nhẹ thôi (Targeted)
            scrambleItems(); 
            
            std::cout << "   -> Fail. Revert & Reduce Rx to " << current_Rx << std::endl;
        }
        iter++;
    }

    items = s_best;
    config.container_size = best_size;
}

void SparrowSolver::compress(double ratio) {
    // Tương tự explore nhưng shrink rate nhỏ dần
    auto start = std::chrono::steady_clock::now();
    std::vector<CompositeShape> s_star = items;
    double best_size = config.container_size;
    
    int loop_iter = 0;
    while (loop_iter < config.max_outer_loops) {
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start).count();
        if (elapsed > config.TLc) break;

        // Dynamic shrink rate
        double progress = (double)elapsed / config.TLc;
        double r = config.Rs_c + (config.Re_c - config.Rs_c) * progress;

        items = s_star; // Luôn bắt đầu từ feasible
        double old_size = config.container_size;
        config.container_size *= (1.0 - r);

        // Reset Weights
        for(auto& row : weights) std::fill(row.begin(), row.end(), 1.0);

        bool feasible = separate(config.Kc, config.Nc, config.TLc - elapsed);

        if (feasible) {
            s_star = items;
            best_size = config.container_size;
             std::cout << "   Compress success -> " << best_size << std::endl;
        } else {
            config.container_size = old_size; // Revert nếu fail
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
    
    // Initial Feasibility check
    int attempts = 0;
    while (!separate(config.Kx, config.Nx * 2, config.TLx) && attempts < 10) {
        std::cout << "[Warn] Initial config infeasible. Expanding..." << std::endl;
        config.container_size *= 1.1;
        initializePlacement(); // Reset vị trí
        attempts++;
    }
    
    if (attempts >= 10) {
        std::cout << "[Error] Could not find initial feasible solution." << std::endl;
        return;
    }

    std::cout << ">>> Initial Feasible Found: " << config.container_size << std::endl;

    double bestSize = config.container_size;
    std::vector<CompositeShape> bestItems = items;
    
    int maxRestarts = 10;
    for (int iter = 0; iter < maxRestarts; ++iter) {
        std::cout << "\n=== Major Iteration " << iter << " ===" << std::endl;
        
        // 1. Descent (Compress)
        descent();
        
        if (config.container_size < bestSize) {
            bestSize = config.container_size;
            bestItems = items;
            std::cout << "   [***] New Record: " << bestSize << std::endl;
        }

        // 2. Explore (Escape local optima)
        // Load lại best để explore từ đó
        items = bestItems;
        config.container_size = bestSize;
        explore();
        
        // Cập nhật lại nếu explore tìm được cái tốt hơn
        if (config.container_size < bestSize) {
            bestSize = config.container_size;
            bestItems = items;
            std::cout << "   [***] Explore found better: " << bestSize << std::endl;
        }
    }

    // Final result
    items = bestItems;
    config.container_size = bestSize;
    std::cout << ">>> Final Size: " << std::setprecision(10) << std::fixed << config.container_size << std::endl;
}