#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "MathUtils.h"

// --- Cấu trúc SoA (Structure of Arrays) ---
// Lưu trữ các hình tròn đại diện (Surrogate Circles)
// Định nghĩa tại đây để CompositeShape có thể chứa nó trực tiếp
struct CirclesSoA {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> r;

    void add(double _x, double _y, double _r) {
        x.push_back(_x);
        y.push_back(_y);
        r.push_back(_r);
    }

    void clear() {
        x.clear();
        y.clear();
        r.clear();
    }

    size_t size() const { return x.size(); }
};

// Hàm tiện ích tính diện tích đa giác (giữ lại từ code cũ)
double getPolygonArea(const std::vector<Vec2>& points);

// --- ConvexPolygon ---
// Đại diện cho một đa giác lồi (thành phần của CompositeShape)
class ConvexPolygon {
private:
    std::vector<Vec2> localVertices; 
    std::vector<Vec2> localNormals;

    // Cache world transform
    mutable std::vector<Vec2> worldVertices; 
    mutable std::vector<Vec2> worldNormals;
    mutable AABB worldAABB;
    
    Transform transform;
    mutable bool isDirty; 

    void updateCache() const;

public:
    ConvexPolygon();
    
    void warmUp() const;
    void setVertices(const std::vector<Vec2>& verts);
    void setTransform(Vec2 pos, double angle);

    const std::vector<Vec2>& getVertices() const;
    const std::vector<Vec2>& getNormals() const;
    const AABB& getAABB() const;
    const std::vector<Vec2>& getLocalVertices() const;

    // Helper: Kiểm tra một điểm có nằm trong đa giác không
    // (Dùng để sinh các hình tròn đại diện)
    bool contains(Vec2 p) const;

    // Metrics
    double diameter = 0.0;
    double convexHullArea = 0.0;
    
    // Hàm tính toán metrics cơ bản (area, diameter)
    void precomputeMetrics();

    bool operator==(const ConvexPolygon& other) const {
        return localVertices == other.localVertices && transform == other.transform;
    }

    bool operator<(const ConvexPolygon& other) const {
        return localVertices < other.localVertices;
    }
};

// --- CompositeShape ---
// Đại diện cho một vật thể phức hợp (chứa nhiều đa giác lồi)
struct CompositeShape {
    std::vector<ConvexPolygon> parts;
    AABB totalAABB;
    Vec2 pos = {0,0};
    double angle = 0.0;

    // --- DATA MỚI CHO THUẬT TOÁN TỐI ƯU ---
    // Lưu trữ các hình tròn đại diện ở hệ tọa độ địa phương (Local Coordinates)
    CirclesSoA circles; 
    
    // ---------------------------------------

    void warmUp() const;
    void addPart(const std::vector<Vec2>& verts);
    void setTransform(Vec2 pos, double angle);

    // Helper: Kiểm tra điểm nằm trong shape
    bool contains(Vec2 p) const;

    // Metrics
    double diameter = 0.0;
    double convexHullArea = 0.0;
    
    // Hàm mới: Sinh các hình tròn đại diện (Voxelization/Medial Axis approximation)
    // quality_scale: 1.0 là mặc định, nhỏ hơn thì sinh ít circle hơn (nhanh hơn nhưng kém chính xác)
    void generateSurrogateCircles(double quality_scale = 1.0);

    bool isComplex() const { return parts.size() > 1; }

    bool operator==(const CompositeShape& other) const {
        if (pos != other.pos || angle != other.angle || parts.size() != other.parts.size()) return false;
        for (size_t i = 0; i < parts.size(); ++i) {
            if (!(parts[i] == other.parts[i])) return false;
        }
        return true;
    }

    bool operator<(const CompositeShape& other) const {
        if (pos != other.pos) return pos < other.pos;
        if (angle != other.angle) return angle < other.angle;
        return parts < other.parts;
    }
};

#endif // GEOMETRY_H