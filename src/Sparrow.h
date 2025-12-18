// Sparrow.h
#ifndef SPARROW_H
#define SPARROW_H

#include <vector>
#include <random>
#include "Geometry.h"

// Cấu trúc lưu trữ cấu hình thuật toán (thêm tham số từ Table 1 trong paper)
struct SparrowConfig {
    double container_size; // Kích thước cạnh hình vuông hiện tại
    int max_iter = 1000;
    int n_samples = 32;    // Số mẫu thử ngẫu nhiên (Algorithm 6)
    int n_threads = -1;

    // THÊM từ Table 1
    int Nc = 100;          // Inner iterations for compress/descent
    int Kc = 5;            // Strikes for compress/descent
    int Nx = 10;           // Inner iterations for explore
    int Kx = 3;            // Strikes for explore
    double TLx = 960.0;    // Time limit explore (seconds)
    double TLc = 240.0;    // Time limit compress (seconds)
    double Rs_c = 0.0005;  // Initial shrink rate compress
    double Re_c = 0.00001; // Final shrink rate compress
    double Rx = 0.001;     // Shrink rate explore
    double Md = 0.95;      // Multiplier down (GLS)
    double Ml = 1.05;      // Multiplier low (GLS)
    double Mu = 1.1;       // Multiplier up (GLS)
    int max_outer_loops = 20;
    int max_restart = 10;
};

class SparrowSolver {
private:
    // --- Data ---
    std::vector<CompositeShape> items;        // Danh sách vật thể hiện tại
    std::vector<std::vector<double>> weights; // Ma trận trọng số GLS (Algorithm 8)
    SparrowConfig config;
    
    std::mt19937 rng; // Bộ sinh số ngẫu nhiên

    // --- Backup State (cho quá trình Descent/Explore) ---
    std::vector<CompositeShape> savedItems;
    double savedContainerSize;

    // --- Internal Helpers ---
    double randomDouble(double min, double max);
    void saveState();
    void loadState();

    // Algorithm 7: Tính chi phí (Overlap + Penalty) tại một vị trí giả định
    double evaluateSample(int itemIdx, Vec2 pos, double angle, const std::vector<CompositeShape>& local_items, const std::vector<std::vector<double>>& local_weights) const;

    // Algorithm 6: Tìm vị trí mới tốt hơn cho một vật thể
    void searchPosition(int itemIdx, std::vector<CompositeShape>& local_items, std::vector<std::vector<double>>& local_weights);

    // Algorithm 8: Cập nhật trọng số Guided Local Search
    void updateWeights(std::vector<std::vector<double>>& local_weights);

    // Added for parallel safety
    void perform_move_items(std::vector<CompositeShape>& local_items, std::vector<std::vector<double>>& local_weights);
    void move_items_multi();
    double calculate_energy(const std::vector<CompositeShape>& local_items) const;

    // New: Better initial placement
    void initializePlacement();

public:
    // Constructor
    SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg);

    // --- Core Algorithms ---

    // Tính tổng năng lượng toàn hệ thống (Overlap vật-vật + Overlap tường)
    double getTotalEnergy();

    // Algorithm 9: Hàm tách rời va chạm (Feasibility Solver)
    // Trả về true nếu tìm được lời giải khả thi (Energy ~ 0)
    bool separate(int kmax, int nmax, double time_limit_sec);  

    // Algorithm 13: Co nhỏ thùng chứa (Scale Down)
    void compress(double ratio);

    void scrambleItems();

    // Algorithm 12: Phá vỡ cấu trúc cục bộ (Perturbation)
    void explore();

    // Algorithm 11: Cố gắng co nhỏ sâu nhất có thể từ trạng thái hiện tại
    void descent();

    // Algorithm 10: Hàm giải chính (Outer Loop)
    void solve();

    // --- Getters ---
    const std::vector<CompositeShape>& getItems() const { return items; }
    double getContainerSize() const { return config.container_size; }
};

#endif // SPARROW_H