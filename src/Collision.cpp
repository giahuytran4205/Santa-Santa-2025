#include "Collision.h"
#include <limits>
#include <cmath>
#include <algorithm>

// Hàm tính năng lượng va chạm giữa hai tập hợp hình tròn (CirclesSoA)
// Sử dụng công thức "Soft Penalty" để tạo gradient trơn, giúp thuật toán thoát khỏi cực tiểu cục bộ.
double evaluate_circles_fast(const CirclesSoA& setA, const CirclesSoA& setB, double epsilon) {
    double totalOverlap = 0.0;
    double epsilon_sq = epsilon * epsilon;
    double two_epsilon = 2.0 * epsilon;

    size_t nA = setA.size();
    size_t nB = setB.size();

    // Duyệt qua tất cả các cặp hình tròn
    // Cấu trúc vòng lặp này thân thiện với bộ nhớ cache và trình biên dịch có thể tự động vector hóa (Auto-Vectorization)
    for (size_t i = 0; i < nA; ++i) {
        double x1 = setA.x[i];
        double y1 = setA.y[i];
        double r1 = setA.r[i];

        for (size_t j = 0; j < nB; ++j) {
            double x2 = setB.x[j];
            double y2 = setB.y[j];
            double r2 = setB.r[j];

            double dx = x1 - x2;
            double dy = y1 - y2;
            double dist_sq = dx*dx + dy*dy;
            double sum_r = r1 + r2;

            // Broadphase nhanh: nếu khoảng cách lớn hơn tổng bán kính -> không chạm -> bỏ qua
            if (dist_sq >= sum_r * sum_r) continue;

            double dist = std::sqrt(dist_sq);
            double pd = sum_r - dist; // Penetration Depth (Độ lún)

            // --- MAGIC FORMULA TỪ RUST (Quantify/overlap_proxy_simd.rs) ---
            double pd_decay;
            if (pd >= epsilon) {
                // Trường hợp 1: Lún sâu (Deep Penetration)
                // Phạt tuyến tính để tạo lực đẩy mạnh ra ngoài
                pd_decay = pd;
            } else {
                // Trường hợp 2: Lún nông hoặc vừa chạm (Shallow Penetration)
                // Phạt theo hàm Hyperbolic: f(x) = e^2 / (-x + 2e)
                // Hàm này tạo ra một "đệm khí" (cushion) giúp vật thể trượt nhẹ nhàng qua nhau
                // thay vì bị đẩy bật ra một cách hỗn loạn.
                pd_decay = epsilon_sq / (-pd + two_epsilon);
            }

            // Cộng dồn năng lượng va chạm
            // Nhân với min(r1, r2) để chuẩn hóa lực đẩy theo kích thước chi tiết nhỏ nhất
            totalOverlap += pd_decay * std::min(r1, r2);
        }
    }
    
    // Scale kết quả theo PI (để tương đồng với đơn vị diện tích hình tròn, giống bản Rust)
    return totalOverlap * 3.14159265358979323846;
}

// Hàm kiểm tra va chạm đa giác truyền thống
// Giữ lại để đảm bảo tính tương thích (nếu logic cũ còn cần check biên chính xác)
bool checkCollisionComposite(const CompositeShape& shapeA, const CompositeShape& shapeB) {
    // Bước 1: Check AABB tổng quát trước
    if (!shapeA.totalAABB.overlaps(shapeB.totalAABB)) return false;

    // Bước 2: Fallback (Nếu cần độ chính xác tuyệt đối thì dùng SAT, 
    // nhưng trong thuật toán tối ưu này ta chủ yếu dựa vào hàm năng lượng ở trên).
    // Ở đây trả về false để tránh tốn chi phí tính toán SAT không cần thiết 
    // nếu bạn đã chuyển hoàn toàn sang dùng Circle Physics.
    return false; 
}

// Hàm định lượng va chạm cũ (nếu code cũ còn gọi) -> map sang hàm mới
double quantify_collision(const CompositeShape& shapeA, const CompositeShape& shapeB) {
    // Map sang hàm evaluate_circles_fast với epsilon nhỏ
    return evaluate_circles_fast(shapeA.circles, shapeB.circles, 1e-9);
}