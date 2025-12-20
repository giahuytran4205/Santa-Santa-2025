// Collision.cpp
#include "Collision.h"
#include <limits>
#include <algorithm>
#include <cmath>

// Helper check SAT (Giữ nguyên logic chuẩn)
bool checkSAT(const ConvexPolygon& polyA, const ConvexPolygon& polyB) {
    // ... (Giữ nguyên code SAT cũ của bạn nếu nó đã đúng logic toán học)
    // Để an toàn, tôi viết lại bản rút gọn chuẩn ở dưới:
    const auto& vertsA = polyA.getVertices();
    const auto& normsA = polyA.getNormals();
    const auto& vertsB = polyB.getVertices();
    const auto& normsB = polyB.getNormals();

    auto isSeparated = [&](const std::vector<Vec2>& norms) {
        for (const auto& axis : norms) {
            double minA = std::numeric_limits<double>::max(), maxA = -minA;
            double minB = std::numeric_limits<double>::max(), maxB = -minB;
            for (const auto& v : vertsA) { double p = v.dot(axis); if(p < minA) minA=p; if(p > maxA) maxA=p; }
            for (const auto& v : vertsB) { double p = v.dot(axis); if(p < minB) minB=p; if(p > maxB) maxB=p; }
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
    double distSq = (shapeA.pos - shapeB.pos).lengthSq();
    double rSum = (shapeA.diameter + shapeB.diameter) * 0.5;
    if (distSq > rSum * rSum) return 0.0;
    
    // 1. Broadphase siêu nhanh
    if (!shapeA.totalAABB.overlaps(shapeB.totalAABB)) return 0.0;

    // QUAN TRỌNG: Giảm epsilon xuống cực nhỏ để các vật có thể xếp khít
    // 1e-5 thay vì 0.01 (1%)
    double R_epsilon = 1e-15;
    double epsilon = R_epsilon * std::max(shapeA.diameter, shapeB.diameter);
    
    // Ngưỡng cắt bỏ: Nếu tách xa quá ngưỡng này thì coi như không va chạm (Energy = 0)
    // Giúp thuật toán không bị nhiễu bởi các vật ở xa nhưng có AABB chồng lấn
    double separation_cutoff = epsilon * 2.0; 

    // Logic tính toán va chạm chi tiết (SAT Projection)
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

            double minOverlap = std::numeric_limits<double>::max();
            bool separated = false;

            // Kiểm tra trục của A
            for (const auto& axis : normsA) {
                double minA = 1e9, maxA = -1e9;
                double minB = 1e9, maxB = -1e9;
                for (auto& v : vertsA) { double p = v.dot(axis); if(p<minA) minA=p; if(p>maxA) maxA=p; }
                for (auto& v : vertsB) { double p = v.dot(axis); if(p<minB) minB=p; if(p>maxB) maxB=p; }
                
                double overlap = std::min(maxA, maxB) - std::max(minA, minB);
                if (overlap < 0) { // Đã tách rời trên trục này
                    separated = true;
                    // Lấy khoảng cách tách rời lớn nhất (ít âm nhất) làm overlap
                    if (overlap > -1e9 && overlap < minOverlap) minOverlap = overlap; 
                    // SAT: Chỉ cần tách rời trên 1 trục là tách rời hoàn toàn
                    // Tuy nhiên ta cần tìm độ tách rời nhỏ nhất (closest distance) nên phải check hết các trục?
                    // Không, với SAT, khoảng cách thực là max(overlap_negative).
                    // Logic cũ của bạn tìm maxNegativeOverlap là đúng.
                } else {
                    if (overlap < minOverlap) minOverlap = overlap;
                }
            }

            // Kiểm tra trục của B (Logic tương tự)
            // ... Để code gọn, ta giả định bạn đã implement đúng việc tìm 
            // 'delta' là khoảng cách xâm nhập (dương) hoặc khoảng cách tách rời (âm)
            // Ở đây tôi viết logic gộp simplified:
            
            // (Đoạn này giả lập lại logic SAT tìm delta chính xác nhất)
            double bestAxisVal = -std::numeric_limits<double>::max(); // Max negative overlap (separation)
            double minPenetration = std::numeric_limits<double>::max(); // Min positive overlap
            bool is_separated_SAT = false;

            auto checkAxes = [&](const std::vector<Vec2>& norms) {
                 for (const auto& axis : norms) {
                    double minA=1e9, maxA=-1e9, minB=1e9, maxB=-1e9;
                    for (auto& v : vertsA) { double p = v.dot(axis); if(p<minA) minA=p; if(p>maxA) maxA=p; }
                    for (auto& v : vertsB) { double p = v.dot(axis); if(p<minB) minB=p; if(p>maxB) maxB=p; }
                    double ov = std::min(maxA, maxB) - std::max(minA, minB);
                    if (ov < 0) {
                        is_separated_SAT = true;
                        if (ov > bestAxisVal) bestAxisVal = ov; // Tìm khoảng hở nhỏ nhất (gần 0 nhất)
                    } else {
                        if (ov < minPenetration) minPenetration = ov;
                    }
                 }
            };
            checkAxes(normsA);
            checkAxes(normsB);

            double delta = is_separated_SAT ? bestAxisVal : minPenetration;

            // --- FIXING THE DECAY LOGIC ---
            
            // Nếu đã tách rời quá xa -> Bỏ qua ngay (Cost = 0)
            if (is_separated_SAT && (-delta > separation_cutoff)) {
                continue; 
            }

            double delta_prime;
            if (delta > epsilon) {
                delta_prime = delta; // Phạt trực tiếp bằng độ lún
            } else {
                // Công thức decay của bài báo
                // Chỉ áp dụng khi gần chạm (-cutoff < delta < epsilon)
                delta_prime = (epsilon * epsilon) / (-delta + 2.0 * epsilon);
            }

            // Cộng dồn vào proxy tổng
            maxProxy += delta_prime * std::min(partA.diameter, partB.diameter);
            interaction_found = true;
        }
    }
    
    if (!interaction_found) return 0.0;

    double lambda_a = std::sqrt(shapeA.convexHullArea);
    double lambda_b = std::sqrt(shapeB.convexHullArea);
    return std::sqrt(maxProxy * std::sqrt(lambda_a * lambda_b));
}