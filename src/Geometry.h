// Geometry.h
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include "MathUtils.h"

double getPolygonArea(const std::vector<Vec2>& points);

// Pole struct
struct Pole {
    Vec2 center;
    double radius;  // Distance to boundary
};

// --- ConvexPolygon ---
class ConvexPolygon {
private:
    std::vector<Vec2> localVertices; 
    std::vector<Vec2> localNormals;

    mutable std::vector<Vec2> worldVertices; 
    mutable std::vector<Vec2> worldNormals;
    mutable AABB worldAABB;
    
    Transform transform;
    mutable bool isDirty; 

    // Hàm nội bộ để cập nhật cache khi cần thiết
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

    // THÊM
    std::vector<Pole> poles;
    double diameter = 0.0;
    double convexHullArea = 0.0;
    void precomputePolesAndMetrics(double precision_rel = 0.01);

    bool operator==(const ConvexPolygon& other) const {
        return localVertices == other.localVertices && transform == other.transform;
    }

    bool operator<(const ConvexPolygon& other) const {
        return localVertices < other.localVertices;  // Lexicographical compare on vertices
    }
};

// --- CompositeShape ---
struct CompositeShape {
    std::vector<ConvexPolygon> parts;
    AABB totalAABB;
    Vec2 pos = {0,0};
    double angle = 0.0;

    void warmUp() const;
    void addPart(const std::vector<Vec2>& verts);
    void setTransform(Vec2 pos, double angle);

    // THÊM
    std::vector<Pole> allPoles;
    double diameter = 0.0;
    double convexHullArea = 0.0;
    void precomputeAllPoles(double precision_rel = 0.01);
    bool isComplex() const { return parts.size() > 1; }  // Hoặc thêm logic total vertices >20 nếu cần

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
        return parts < other.parts;  // Lexicographical on parts
    }
};

#endif // GEOMETRY_H