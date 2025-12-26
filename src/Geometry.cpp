#include "Geometry.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <queue>
#include <vector>

// ==========================================
// HELPER FUNCTIONS & MATH UTILS
// ==========================================

double getPolygonArea(const std::vector<Vec2>& points) {
    double area = 0.0;
    size_t n = points.size();
    if (n < 3) return 0.0;
    for (size_t i = 0; i < n; ++i) {
        area += (points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y);
    }
    return std::abs(area) * 0.5;
}

// Tính khoảng cách từ điểm P đến đoạn thẳng AB
// Trả về khoảng cách bình phương để tối ưu
double distToSegmentSq(Vec2 p, Vec2 a, Vec2 b) {
    Vec2 ab = b - a;
    double lenSq = ab.lengthSq();
    if (lenSq < 1e-12) return p.distSq(a);
    double t = std::max(0.0, std::min(1.0, (p - a).dot(ab) / lenSq));
    return p.distSq(a + ab * t);
}

// Lấy điểm gần nhất trên đoạn thẳng AB tính từ P
Vec2 getClosestPointOnSegment(Vec2 p, Vec2 a, Vec2 b) {
    Vec2 ab = b - a;
    double lenSq = ab.lengthSq();
    if (lenSq < 1e-12) return a;
    double t = std::max(0.0, std::min(1.0, (p - a).dot(ab) / lenSq));
    return a + ab * t;
}

// ==========================================
// CONVEX POLYGON
// ==========================================

ConvexPolygon::ConvexPolygon() : isDirty(true), diameter(0), convexHullArea(0) {}

void ConvexPolygon::setVertices(const std::vector<Vec2>& verts) {
    localVertices = verts;
    isDirty = true;
    precomputeMetrics();
}

void ConvexPolygon::setTransform(Vec2 pos, double angle) {
    if (transform.pos != pos || transform.angle != angle) {
        transform.pos = pos;
        transform.angle = angle;
        isDirty = true;
    }
}

void ConvexPolygon::updateCache() const {
    if (!isDirty) return;
    worldVertices.resize(localVertices.size());
    
    double c = std::cos(transform.angle);
    double s = std::sin(transform.angle);
    double minX = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();

    for (size_t i = 0; i < localVertices.size(); ++i) {
        double x = localVertices[i].x * c - localVertices[i].y * s + transform.pos.x;
        double y = localVertices[i].x * s + localVertices[i].y * c + transform.pos.y;
        worldVertices[i] = {x, y};
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }
    worldAABB = { {minX, minY}, {maxX, maxY} };
    isDirty = false;
}

void ConvexPolygon::warmUp() const { updateCache(); }
const std::vector<Vec2>& ConvexPolygon::getVertices() const { updateCache(); return worldVertices; }
const std::vector<Vec2>& ConvexPolygon::getNormals() const { return worldNormals; } 
const AABB& ConvexPolygon::getAABB() const { updateCache(); return worldAABB; }
const std::vector<Vec2>& ConvexPolygon::getLocalVertices() const { return localVertices; }

// Ray Casting check point inside
bool ConvexPolygon::contains(Vec2 p) const {
    const auto& verts = getVertices();
    bool inside = false;
    size_t n = verts.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        if (((verts[i].y > p.y) != (verts[j].y > p.y)) &&
            (p.x < (verts[j].x - verts[i].x) * (p.y - verts[i].y) / (verts[j].y - verts[i].y) + verts[i].x)) {
            inside = !inside;
        }
    }
    return inside;
}

void ConvexPolygon::precomputeMetrics() {
    convexHullArea = getPolygonArea(localVertices);
    diameter = 0.0;
    for (size_t i = 0; i < localVertices.size(); ++i) {
        for (size_t j = i + 1; j < localVertices.size(); ++j) {
            double d = localVertices[i].distSq(localVertices[j]);
            if (d > diameter) diameter = d;
        }
    }
    diameter = std::sqrt(diameter);
}

// ==========================================
// COMPOSITE SHAPE
// ==========================================

void CompositeShape::addPart(const std::vector<Vec2>& verts) {
    ConvexPolygon poly;
    poly.setVertices(verts);
    parts.push_back(poly);
    convexHullArea += poly.convexHullArea; // Approximate
    if (poly.diameter > diameter) diameter = poly.diameter; 
}

void CompositeShape::setTransform(Vec2 p, double a) {
    if (pos != p || angle != a) {
        pos = p;
        angle = a;
        double minX = std::numeric_limits<double>::max();
        double maxX = -std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        double maxY = -std::numeric_limits<double>::max();
        for (auto& part : parts) {
            part.setTransform(pos, angle);
            const AABB& b = part.getAABB();
            if (b.min.x < minX) minX = b.min.x;
            if (b.max.x > maxX) maxX = b.max.x;
            if (b.min.y < minY) minY = b.min.y;
            if (b.max.y > maxY) maxY = b.max.y;
        }
        totalAABB = { {minX, minY}, {maxX, maxY} };
    }
}

void CompositeShape::warmUp() const { for (const auto& part : parts) part.warmUp(); }

bool CompositeShape::contains(Vec2 p) const {
    if (!totalAABB.contains(p)) return false;
    for (const auto& part : parts) {
        if (part.getAABB().contains(p)) {
            if (part.contains(p)) return true;
        }
    }
    return false;
}

// ==========================================
// THUẬT TOÁN UNIFIED POLE OF INACCESSIBILITY
// ==========================================

// Cấu trúc Node cho Quadtree search (giống POINode của Rust)
struct POINode {
    Vec2 center;
    double h; // Half-size of the cell (width/2)
    double dist; // Fitness value (distance to boundary/poles)

    // Max potential distance for any point in this cell
    double getUpperBound() const {
        return dist + h * 1.41421356; // dist + h * sqrt(2)
    }

    bool operator<(const POINode& other) const {
        return dist < other.dist; // Max-heap based on distance
    }
};

struct CircleData { double x, y, r; };

// Kiểm tra xem điểm Q có bị "che" bởi phần nào khác của vật thể không
// Dùng để loại bỏ các cạnh nội bộ (internal edges)
bool isPointInternal(const CompositeShape& shape, Vec2 q, size_t ignorePartIdx) {
    // Check điểm q có nằm TRONG bất kỳ part nào khác không (không tính biên)
    // Dung sai nhỏ để tránh nhiễu số học
    const double EPS = 1e-4; 
    for (size_t i = 0; i < shape.parts.size(); ++i) {
        if (i == ignorePartIdx) continue;
        const auto& part = shape.parts[i];
        
        // Nếu part này chứa q, nghĩa là q nằm bên trong hợp của shape
        // -> Cạnh sinh ra q là cạnh nội bộ.
        // Dùng ray casting check
        if (part.getAABB().contains(q) && part.contains(q)) {
            // Check kỹ hơn khoảng cách biên để chắc chắn nó nằm sâu bên trong
            // hoặc chồng lấn thực sự
            return true;
        }
    }
    return false;
}

// Tính khoảng cách từ P đến "Biên bao" (Unified Boundary) của CompositeShape
// Logic: Khoảng cách ngắn nhất tới bất kỳ cạnh nào mà cạnh đó KHÔNG phải là cạnh nội bộ.
double getDistToUnifiedBoundary(const CompositeShape& shape, Vec2 p) {
    // 1. Kiểm tra P có nằm trong Shape không (bất kỳ phần nào)
    bool inside = false;
    for (const auto& part : shape.parts) {
        if (part.getAABB().contains(p) && part.contains(p)) {
            inside = true;
            break;
        }
    }
    // Nếu nằm ngoài -> trả về số âm (hoặc -1)
    if (!inside) return -1.0;

    double minD = std::numeric_limits<double>::max();

    // 2. Tìm khoảng cách tới cạnh gần nhất
    for (size_t i = 0; i < shape.parts.size(); ++i) {
        const auto& part = shape.parts[i];
        const auto& verts = part.getVertices();
        size_t n = verts.size();

        for (size_t k = 0; k < n; ++k) {
            Vec2 v1 = verts[k];
            Vec2 v2 = verts[(k + 1) % n];
            
            // Tìm điểm gần nhất Q trên đoạn thẳng v1-v2
            Vec2 q = getClosestPointOnSegment(p, v1, v2);
            
            // QUAN TRỌNG: Kiểm tra Q có phải là điểm "nội bộ" (bị che bởi part khác) không?
            // Nếu Q nằm trong một part khác -> Cạnh này là đường nối bên trong -> Bỏ qua
            if (isPointInternal(shape, q, i)) {
                continue; 
            }

            double dSq = p.distSq(q);
            if (dSq < minD) minD = dSq;
        }
    }
    
    return std::sqrt(minD);
}

// Tìm Pole tiếp theo dựa trên Quadtree
CircleData computePole(const CompositeShape& shape, const std::vector<CircleData>& existingPoles) {
    const double precision = 1e-3; // Độ chính xác dừng
    AABB rootBox = shape.totalAABB;
    
    // Hàm đánh giá Fitness cho điểm P
    auto evaluate = [&](Vec2 p) -> double {
        // Dist to unified boundary (Dương nếu trong, Âm nếu ngoài)
        double dPoly = getDistToUnifiedBoundary(shape, p);
        if (dPoly <= 0) return dPoly; // Nằm ngoài -> Invalid

        // Dist to existing poles (signed distance)
        // Chúng ta muốn điểm P nằm NGOÀI các pole đã có
        // Dist = Khoảng cách tới tâm Pole - Bán kính Pole
        double dPoles = std::numeric_limits<double>::max();
        for (const auto& pole : existingPoles) {
            // Signed distance function của hình tròn:
            // d > 0: Nằm ngoài pole
            // d < 0: Nằm trong pole (phạt)
            double d = p.dist({pole.x, pole.y}) - pole.r;
            if (d < dPoles) dPoles = d;
        }

        // Fitness là min của (khoảng cách tới biên, khoảng cách tới các pole khác)
        // Điều này ép hình tròn mới phải vừa nằm trong shape, vừa không đè lên pole cũ
        return std::min(dPoly, dPoles);
    };

    std::priority_queue<POINode> queue;
    
    // Init Root Node
    double cellSize = std::max(rootBox.width(), rootBox.height()) * 0.5;
    Vec2 center = rootBox.center();
    queue.push({center, cellSize, evaluate(center)});

    CircleData bestPole = {center.x, center.y, 0.0}; // Default bad pole

    while (!queue.empty()) {
        POINode node = queue.top();
        queue.pop();

        // Pruning: Nếu node này không thể chứa điểm tốt hơn best hiện tại
        if (node.getUpperBound() <= bestPole.r) continue;

        // Cập nhật Best
        if (node.dist > bestPole.r) {
            bestPole = {node.center.x, node.center.y, node.dist};
        }

        // Split
        if (node.h > precision) {
            double h = node.h * 0.5;
            // 4 children centered at offsets
            double offsets[4][2] = {{-1, -1}, {1, -1}, {-1, 1}, {1, 1}};
            for (int i = 0; i < 4; ++i) {
                Vec2 c = {node.center.x + offsets[i][0] * h, node.center.y + offsets[i][1] * h};
                // Quick check: Nếu bbox con nằm hoàn toàn ngoài AABB shape -> Bỏ qua (Optimization)
                // (Ở đây ta check đơn giản bằng evaluate tại tâm)
                queue.push({c, h, evaluate(c)});
            }
        }
    }

    return bestPole;
}

// --- MAIN FUNCTION: GENERATE SURROGATE CIRCLES ---
// Logic: Iterative Pole of Inaccessibility
void CompositeShape::generateSurrogateCircles(double quality_scale) {
    circles.clear();
    
    // 1. Reset transform về local (0,0)
    Vec2 oldPos = pos;
    double oldAngle = angle;
    setTransform({0, 0}, 0);
    
    std::vector<CircleData> poles;
    
    // Số lượng Pole mục tiêu (giống Rust config)
    // quality_scale = 1.0 -> khoảng 16-32 poles là đẹp
    int maxPoles = (int)(24 * quality_scale);
    if (maxPoles < 4) maxPoles = 4;

    // Diện tích mục tiêu phủ (Heuristic)
    double targetArea = 0.0;
    for(const auto& p : parts) targetArea += p.convexHullArea;
    targetArea *= 0.95; // Phủ 95% diện tích

    double currentArea = 0.0;

    for (int i = 0; i < maxPoles; ++i) {
        // Tìm Pole tốt nhất
        CircleData p = computePole(*this, poles);
        
        // Nếu bán kính quá nhỏ (vào khe hẹp quá mức) hoặc âm (không tìm thấy chỗ trống) -> Dừng
        if (p.r < 1e-3) break;
        
        poles.push_back(p);
        circles.add(p.x, p.y, p.r);
        
        currentArea += 3.14159 * p.r * p.r;
        // Break sớm nếu đã phủ kín (tương đối)
        if (currentArea > targetArea) break;
    }
    
    // FALLBACK: Nếu shape quá dị, ít nhất có 1 circle ở tâm
    if (circles.size() == 0) {
        Vec2 c = totalAABB.center();
        double r = std::min(totalAABB.width(), totalAABB.height()) * 0.25;
        circles.add(c.x, c.y, r);
    }

    // 2. Restore transform
    setTransform(oldPos, oldAngle);
}