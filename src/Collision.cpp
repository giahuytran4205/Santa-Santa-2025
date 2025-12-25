// Collision.cpp
#include "Collision.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include <immintrin.h> // Thư viện SIMD (AVX/AVX2)

// --- Helper: Tính toán Min/Max Projection bằng AVX2 ---
// Hàm này thay thế cho vòng lặp for duyệt đỉnh thông thường
// Giả định Vec2 gồm {double x, double y} nằm liền kề trong bộ nhớ
void project_AVX2(const std::vector<Vec2>& verts, Vec2 axis, double& minOut, double& maxOut) {
    size_t n = verts.size();
    const double* data = reinterpret_cast<const double*>(verts.data());

    double minVal = std::numeric_limits<double>::max();
    double maxVal = -std::numeric_limits<double>::max();

    size_t i = 0;

    // Kiểm tra hỗ trợ AVX2 tại thời điểm biên dịch (thường được bật bởi -march=native)
#if defined(__AVX2__)
    // 1. Chuẩn bị thanh ghi chứa Axis: [nx, ny, nx, ny]
    // Vì mỗi Vec2 là (x, y), ta cần nhân x*nx và y*ny.
    __m256d axisVec = _mm256_setr_pd(axis.x, axis.y, axis.x, axis.y);
    
    // Khởi tạo min/max vectors với giá trị vô cực
    __m256d vMin = _mm256_set1_pd(minVal);
    __m256d vMax = _mm256_set1_pd(maxVal);

    // 2. Vòng lặp chính: Xử lý 4 đỉnh mỗi lần (2 thanh ghi AVX)
    // Mỗi thanh ghi AVX2 256-bit chứa được 4 số double -> Tương đương 2 struct Vec2.
    // Ta xử lý 2 thanh ghi một lúc để tận dụng pipeline của CPU (Loop Unrolling).
    for (; i + 3 < n; i += 4) {
        // Tải 2 đỉnh đầu tiên (x0, y0, x1, y1)
        __m256d v0 = _mm256_loadu_pd(data + i * 2); 
        // Tải 2 đỉnh tiếp theo (x2, y2, x3, y3)
        __m256d v1 = _mm256_loadu_pd(data + (i + 2) * 2);

        // Nhân với axis: [x0*nx, y0*ny, x1*nx, y1*ny]
        __m256d mul0 = _mm256_mul_pd(v0, axisVec);
        __m256d mul1 = _mm256_mul_pd(v1, axisVec);

        // Cộng ngang (Horizontal Add) để ra Dot Product
        // _mm256_hadd_pd(A, A) với A=[p0x, p0y, p1x, p1y] -> [p0, p1, p0, p1]
        // Kết quả dot0 chứa 2 giá trị dot product (lặp lại)
        __m256d dot0 = _mm256_hadd_pd(mul0, mul0);
        __m256d dot1 = _mm256_hadd_pd(mul1, mul1);

        // Cập nhật Min/Max song song
        vMin = _mm256_min_pd(vMin, dot0);
        vMin = _mm256_min_pd(vMin, dot1);
        vMax = _mm256_max_pd(vMax, dot0);
        vMax = _mm256_max_pd(vMax, dot1);
    }

    // 3. Gom kết quả từ thanh ghi AVX về scalar
    // vMin đang chứa [min01, min01, min23, min23]
    double tempMin[4], tempMax[4];
    _mm256_storeu_pd(tempMin, vMin);
    _mm256_storeu_pd(tempMax, vMax);

    for(int k=0; k<4; ++k) {
        if(tempMin[k] < minVal) minVal = tempMin[k];
        if(tempMax[k] > maxVal) maxVal = tempMax[k];
    }
#endif

    // 4. Xử lý phần dư (Tail loop) bằng Scalar
    // Xử lý các đỉnh còn lại nếu n không chia hết cho 4 (hoặc nếu không có AVX2)
    for (; i < n; ++i) {
        double p = verts[i].dot(axis);
        if (p < minVal) minVal = p;
        if (p > maxVal) maxVal = p;
    }

    minOut = minVal;
    maxOut = maxVal;
}

// Helper check SAT (Đã tối ưu SIMD)
bool checkSAT(const ConvexPolygon& polyA, const ConvexPolygon& polyB) {
    const auto& vertsA = polyA.getVertices();
    const auto& normsA = polyA.getNormals();
    const auto& vertsB = polyB.getVertices();
    const auto& normsB = polyB.getNormals();

    auto isSeparated = [&](const std::vector<Vec2>& norms) {
        for (const auto& axis : norms) {
            double minA, maxA, minB, maxB;
            // Sử dụng SIMD Projection
            project_AVX2(vertsA, axis, minA, maxA);
            project_AVX2(vertsB, axis, minB, maxB);

            if (maxA < minB || maxB < minA) return true; // Separated
        }
        return false;
    };

    if (isSeparated(normsA)) return false;
    if (isSeparated(normsB)) return false;
    return true;
}

bool checkCollisionComposite(const CompositeShape& shapeA, const CompositeShape& shapeB) {
    if (!shapeA.totalAABB.overlaps(shapeB.totalAABB)) return false;
    for (const auto& partA : shapeA.parts) {
        if (!partA.getAABB().overlaps(shapeB.totalAABB)) continue;
        for (const auto& partB : shapeB.parts) {
            if (!partA.getAABB().overlaps(partB.getAABB())) continue;
            if (checkSAT(partA, partB)) return true;
        }
    }
    return false;
}

double quantify_collision(const CompositeShape& shapeA, const CompositeShape& shapeB) {
    // 1. Broadphase siêu nhanh
    if (!shapeA.totalAABB.overlaps(shapeB.totalAABB)) return 0.0;

    // QUAN TRỌNG: Giảm epsilon xuống cực nhỏ để các vật có thể xếp khít
    // 1e-15 như yêu cầu của bạn (cần độ chính xác double)
    double R_epsilon = 1e-15;
    double epsilon = R_epsilon * std::max(shapeA.diameter, shapeB.diameter);
    
    // Ngưỡng cắt bỏ
    double separation_cutoff = epsilon * 2.0; 

    double maxProxy = 0.0;
    bool interaction_found = false;

    for (const auto& partA : shapeA.parts) {
        for (const auto& partB : shapeB.parts) {
            // Check AABB con
            if (!partA.getAABB().overlaps(partB.getAABB())) continue;

            const auto& vertsA = partA.getVertices();
            const auto& normsA = partA.getNormals();
            const auto& vertsB = partB.getVertices();
            const auto& normsB = partB.getNormals();

            double bestAxisVal = -std::numeric_limits<double>::max(); // Max negative overlap (separation)
            double minPenetration = std::numeric_limits<double>::max(); // Min positive overlap
            bool is_separated_SAT = false;

            // Lambda kiểm tra trục (đã tích hợp SIMD)
            auto checkAxes = [&](const std::vector<Vec2>& norms) {
                 for (const auto& axis : norms) {
                    double minA, maxA, minB, maxB;
                    
                    // --- SIMD ACCELERATION HERE ---
                    project_AVX2(vertsA, axis, minA, maxA);
                    project_AVX2(vertsB, axis, minB, maxB);
                    // ------------------------------

                    double ov = std::min(maxA, maxB) - std::max(minA, minB);
                    
                    if (ov < 0) {
                        is_separated_SAT = true;
                        // Tìm khoảng tách rời nhỏ nhất (gần 0 nhất, tức là ov lớn nhất trong số các số âm)
                        if (ov > bestAxisVal) bestAxisVal = ov; 
                    } else {
                        if (ov < minPenetration) minPenetration = ov;
                    }
                 }
            };

            checkAxes(normsA);
            checkAxes(normsB);

            double delta = is_separated_SAT ? bestAxisVal : minPenetration;

            // --- DECAY LOGIC (Giữ nguyên) ---
            
            // Nếu đã tách rời quá xa -> Bỏ qua ngay
            if (is_separated_SAT && (-delta > separation_cutoff)) {
                continue; 
            }

            double delta_prime;
            if (delta > epsilon) {
                delta_prime = delta; // Phạt trực tiếp bằng độ lún
            } else {
                // Công thức decay (Hyperbolic repulsion)
                delta_prime = (epsilon * epsilon) / (-delta + 2.0 * epsilon);
            }

            maxProxy += delta_prime * std::min(partA.diameter, partB.diameter);
            interaction_found = true;
        }
    }
    
    if (!interaction_found) return 0.0;

    double lambda_a = std::sqrt(shapeA.convexHullArea);
    double lambda_b = std::sqrt(shapeB.convexHullArea);
    return std::sqrt(maxProxy * std::sqrt(lambda_a * lambda_b));
}