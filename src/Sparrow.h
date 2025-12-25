#ifndef SPARROW_H
#define SPARROW_H

#include <vector>
#include <random>
#include <chrono>
#include "Geometry.h"

// Cấu trúc cấu hình (Giữ nguyên các tham số tối ưu)
struct SparrowConfig {
    double container_size;
    int max_iter = 1000;
    int n_samples = 1024;  // Tăng mẫu thử vì worker chạy rất nhanh
    int n_threads = -1;

    // Tham số thuật toán
    int Nc = 100;          
    int Kc = 5;            
    int Nx = 20;           
    int Kx = 3;            
    double TLx = 60.0;     // Thời gian explore
    double TLc = 60.0;     // Thời gian compress
    double Rs_c = 0.002;
    double Re_c = 0.0001;
    double Rx = 0.005;     // Bước nhảy nhỏ để chính xác
    double Md = 0.95;
    double Ml = 1.05;
    double Mu = 1.2;
    int max_outer_loops = 50;
    int max_restart = 10;
    
    void loadFromFile(const std::string& filename); 
};

// --- WORKER CLASS ---
// Mô phỏng struct SeparatorWorker trong Rust
// Mỗi worker sở hữu bộ nhớ riêng để tránh race condition và false sharing
class SeparatorWorker {
public:
    int id;
    std::vector<CompositeShape> items;        // Local state (Bản sao cục bộ)
    std::vector<std::vector<double>> weights; // Local weights (Bản sao cục bộ)
    std::mt19937 rng;                         // Local RNG (Độc lập hoàn toàn)
    double current_energy;

    // Buffer tái sử dụng để tránh cấp phát động liên tục
    std::vector<std::pair<double, Transform>> candidate_buffer; 

    // Constructor
    SeparatorWorker(int id, int seed, const std::vector<CompositeShape>& initItems, size_t n);

    // Đồng bộ trạng thái từ Master (Rust: worker.load)
    void load(const std::vector<CompositeShape>& masterItems, const std::vector<std::vector<double>>& masterWeights);

    // Thực hiện di chuyển các vật thể va chạm (Rust: worker.move_items)
    void move_items(const SparrowConfig& config);

private:
    // Các hàm nội bộ của worker (Lock-free)
    void searchPosition(int itemIdx, const SparrowConfig& config);
    double evaluateSample(int itemIdx, Vec2 pos, double angle, const SparrowConfig& config);
    double calculate_total_energy(const SparrowConfig& config);
    double randomDouble(double min, double max);
};

// --- MASTER SOLVER ---
class SparrowSolver {
private:
    std::vector<CompositeShape> items;        // Master State (Global Truth)
    std::vector<std::vector<double>> weights; // Master Weights
    SparrowConfig config;
    std::mt19937 rng;

    // Danh sách các công nhân kiên định
    std::vector<SeparatorWorker> workers; 

    // Backup State
    std::vector<CompositeShape> savedItems;
    double savedContainerSize;

    void saveState();
    void loadState();
    void initializePlacement();
    
    // Hàm cập nhật trọng số GLS trên Master (Algorithm 8)
    void updateMasterWeights();
    double getTotalEnergy() const;

public:
    SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg);

    // Algorithm 9: Tách rời va chạm (Parallel Fork-Join)
    bool separate(int kmax, int nmax, std::chrono::steady_clock::time_point deadline);  
    
    void compress(double ratio);
    void scrambleItems();
    void explore();
    void descent();
    void solve();

    const std::vector<CompositeShape>& getItems() const { return items; }
    double getContainerSize() const { return config.container_size; }
};

#endif // SPARROW_H