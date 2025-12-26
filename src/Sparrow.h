#ifndef SPARROW_H
#define SPARROW_H

#include <vector>
#include <random>
#include <chrono>
#include <string>
#include "Geometry.h"
#include "Collision.h" 

// Cấu trúc cấu hình
struct SparrowConfig {
    double container_size = 1000.0;
    int max_iter = 1000;
    int n_samples = 2048;  // Tăng số lượng mẫu thử (samples) vì thuật toán check va chạm mới rất nhanh
    int n_threads = -1;

    // Tham số thuật toán (giữ nguyên các tham số tối ưu)
    int Nc = 100;
    int Kc = 5;
    int Nx = 20;
    int Kx = 3;
    double TLx = 60.0;     // Thời gian explore
    double TLc = 60.0;     // Thời gian compress
    double Rs_c = 0.002;
    double Re_c = 0.0001;
    double Rx = 0.005;
    double Md = 0.95;
    double Ml = 1.05;
    double Mu = 1.2;
    int max_outer_loops = 50;
    int max_restart = 10;

    void loadFromFile(const std::string& filename);
};

// --- WORKER CLASS ---
// Mỗi worker chạy trên một luồng riêng để đánh giá và di chuyển vật thể
// Sử dụng bộ nhớ cục bộ (Local State) để tránh Race Condition và False Sharing
class SeparatorWorker {
public:
    int id;
    std::vector<CompositeShape> items;        // Local state (Bản sao cục bộ của các vật thể)
    std::vector<std::vector<double>> weights; // Local weights (Trọng số phạt va chạm)
    std::mt19937 rng;                         // Local RNG (Sinh số ngẫu nhiên độc lập)
    double current_energy;

    // Buffer tái sử dụng để lưu các vị trí ứng viên trong quá trình tìm kiếm
    // Giúp tránh cấp phát bộ nhớ động liên tục (tối ưu hiệu năng)
    std::vector<std::pair<double, Transform>> candidate_buffer;

    SeparatorWorker(int id, int seed, const std::vector<CompositeShape>& initItems, size_t n);

    // Đồng bộ trạng thái từ Master xuống Worker
    void load(const std::vector<CompositeShape>& masterItems, const std::vector<std::vector<double>>& masterWeights);
    
    // Thực hiện vòng lặp di chuyển các vật thể va chạm
    void move_items(const SparrowConfig& config);

private:
    // Tìm kiếm vị trí mới tốt hơn cho một vật thể
    void searchPosition(int itemIdx, const SparrowConfig& config);
    
    // Đánh giá năng lượng (mức độ va chạm) tại một vị trí cụ thể
    // Hàm này sẽ sử dụng logic Circle Physics mới
    double evaluateSample(int itemIdx, Vec2 pos, double angle, const SparrowConfig& config);
    
    // Tính tổng năng lượng của toàn bộ hệ thống (dùng để so sánh kết quả giữa các worker)
    double calculate_total_energy(const SparrowConfig& config);
    
    double randomDouble(double min, double max);
};

// --- MASTER SOLVER ---
// Quản lý trạng thái toàn cục và điều phối các Worker
class SparrowSolver {
private:
    std::vector<CompositeShape> items;        // Master State (Sự thật duy nhất)
    std::vector<std::vector<double>> weights; // Master Weights (Trọng số GLS)
    SparrowConfig config;
    std::mt19937 rng;

    std::vector<SeparatorWorker> workers;

    // Backup State để hồi phục khi thuật toán đi vào ngõ cụt
    std::vector<CompositeShape> savedItems;
    double savedContainerSize;

    void saveState();
    void loadState();
    
    // --- THAY ĐỔI CHÍNH: CHIẾN LƯỢC KHỞI TẠO LBF ---
    // Thay thế hàm initializePlacement() ngẫu nhiên cũ bằng thuật toán xây dựng thông minh
    // LBF (Least Bad Fit): Xếp từng vật vào vị trí ít tệ nhất ngay từ đầu.
    void constructLBF();
    
    // Hàm sinh hình tròn đại diện (Surrogate Circles) cho toàn bộ vật thể
    // Được gọi một lần duy nhất khi khởi tạo Solver
    void generateSurrogateCircles();

    // Cập nhật trọng số phạt dựa trên va chạm (Guided Local Search)
    void updateMasterWeights();
    
    double getTotalEnergy() const;

public:
    SparrowSolver(const std::vector<CompositeShape>& inputItems, SparrowConfig cfg);

    // Thuật toán tách rời (Algorithm 9) chạy song song
    bool separate(int kmax, int nmax, std::chrono::steady_clock::time_point deadline);

    // Các chiến lược tối ưu hóa cấp cao
    void compress(double ratio);
    void scrambleItems(); // Xáo trộn khi bị kẹt
    void explore();       // Tìm kiếm không gian mới
    void descent();       // Giảm dần kích thước
    
    // Hàm chính để chạy toàn bộ quy trình
    void solve();

    const std::vector<CompositeShape>& getItems() const { return items; }
    double getContainerSize() const { return config.container_size; }
};

#endif // SPARROW_H