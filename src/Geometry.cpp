// Geometry.cpp
#include "Geometry.h"
#include <limits>
#include <queue>
#include <functional>  // Cho CellCompare

// Helper: Tính diện tích đa giác (Shoelace formula) - Copy từ SVGExport để tự chứa
double getPolygonArea(const std::vector<Vec2>& points) {
    double area = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        size_t j = (i + 1) % points.size();
        area += points[i].x * points[j].y;
        area -= points[j].x * points[i].y;
    }
    return std::abs(area) / 2.0;
}

// Helpers cho poles
double distToLine(Vec2 p, Vec2 a, Vec2 b) {  // Distance point to line segment
    Vec2 ab = b - a;
    double len2 = ab.dot(ab);
    if (len2 == 0) return (p - a).length();
    double t = std::max(0.0, std::min(1.0, (p - a).dot(ab) / len2));
    Vec2 proj = a + ab * t;
    return (p - proj).length();
}

bool isPointInsidePolygon(Vec2 p, const std::vector<Vec2>& verts) {
    int n = verts.size();
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        Vec2 v1 = verts[i] - p;
        Vec2 v2 = verts[(i+1)%n] - p;
        sum += std::atan2(v1.x * v2.y - v1.y * v2.x, v1.dot(v2));
    }
    return std::abs(sum) > PI;  // Winding >0 for inside
}

double signedDistToPolygon(Vec2 p, const std::vector<Vec2>& verts) {
    double minDist = std::numeric_limits<double>::max();
    int n = verts.size();
    for (int i = 0; i < n; ++i) {
        minDist = std::min(minDist, distToLine(p, verts[i], verts[(i+1)%n]));
    }
    bool inside = isPointInsidePolygon(p, verts);
    return inside ? minDist : -minDist;
}

// Quadtree Cell
struct Cell {
    Vec2 center;
    double halfSize;
    double dist;  // signed dist of center
    double potential() const { return dist + halfSize * std::sqrt(2.0); }  // max possible dist in cell
};

struct CellCompare {
    bool operator()(const Cell& a, const Cell& b) { return a.potential() < b.potential(); }  // Max heap
};

// --- ConvexPolygon Implementation ---

ConvexPolygon::ConvexPolygon() : isDirty(true) {}

void ConvexPolygon::warmUp() const { updateCache(); }

void ConvexPolygon::updateCache() const {
    if (!isDirty) return;

    double minX = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();

    size_t count = localVertices.size();
    
    // 1. Transform Vertices & Compute AABB
    for (size_t i = 0; i < count; ++i) {
        Vec2 v = localVertices[i].rotate(transform.sin_a, transform.cos_a) + transform.pos;
        worldVertices[i] = v;

        if (v.x < minX) minX = v.x;
        if (v.x > maxX) maxX = v.x;
        if (v.y < minY) minY = v.y;
        if (v.y > maxY) maxY = v.y;
    }
    worldAABB = {{minX, minY}, {maxX, maxY}};

    // 2. Transform Normals
    for (size_t i = 0; i < count; ++i) {
        worldNormals[i] = localNormals[i].rotate(transform.sin_a, transform.cos_a);
    }

    isDirty = false;
}

void ConvexPolygon::setVertices(const std::vector<Vec2>& verts) {
    localVertices = verts;
    size_t count = verts.size();
    
    worldVertices.resize(count);
    worldNormals.resize(count);
    localNormals.resize(count);

    // Pre-compute Local Normals
    for (size_t i = 0; i < count; ++i) {
        Vec2 p1 = localVertices[i];
        Vec2 p2 = localVertices[(i + 1) % count];
        Vec2 edge = p2 - p1;
        // Pháp tuyến (-y, x)
        Vec2 normal = {-edge.y, edge.x};
        double len = std::sqrt(normal.x*normal.x + normal.y*normal.y);
        localNormals[i] = normal * (1.0 / len);
    }
    
    isDirty = true;
}

void ConvexPolygon::setTransform(Vec2 pos, double angle) {
    transform.set(pos, angle);
    isDirty = true; 
}

const std::vector<Vec2>& ConvexPolygon::getVertices() const {
    updateCache();
    return worldVertices;
}

const std::vector<Vec2>& ConvexPolygon::getNormals() const {
    updateCache();
    return worldNormals;
}

const AABB& ConvexPolygon::getAABB() const {
    updateCache();
    return worldAABB;
}

const std::vector<Vec2>& ConvexPolygon::getLocalVertices() const {
    return localVertices;
}

void ConvexPolygon::precomputePolesAndMetrics(double precision_rel) {
    convexHullArea = getPolygonArea(localVertices);
    diameter = 0.0;
    for (size_t i = 0; i < localVertices.size(); ++i) {
        for (size_t j = i + 1; j < localVertices.size(); ++j) {
            double d = (localVertices[i] - localVertices[j]).length();
            diameter = std::max(diameter, d);
        }
    }
    double precision = precision_rel * diameter;

    // Polylabel quadtree for multiple poles
    poles.clear();
    double minX = std::numeric_limits<double>::max(), maxX = -minX;
    double minY = std::numeric_limits<double>::max(), maxY = -minY;
    for (const auto& v : localVertices) {
        minX = std::min(minX, v.x); maxX = std::max(maxX, v.x);
        minY = std::min(minY, v.y); maxY = std::max(maxY, v.y);
    }
    double cx = (minX + maxX) / 2.0;
    double cy = (minY + maxY) / 2.0;
    double halfSize = std::max(maxX - minX, maxY - minY) / 2.0;
    Cell initial = { {cx, cy}, halfSize, signedDistToPolygon({cx, cy}, localVertices) };

    std::priority_queue<Cell, std::vector<Cell>, CellCompare> pq;
    pq.push(initial);

    double bestDist = 0.0;
    while (!pq.empty()) {
        Cell cell = pq.top(); pq.pop();
        if (cell.dist < 0) continue;  // Outside

        if (cell.halfSize < precision) {
            poles.push_back({cell.center, cell.dist});
            bestDist = std::max(bestDist, cell.dist);
            continue;
        }

        if (cell.potential() <= bestDist) continue;  // Prune

        // Subdivide
        double h = cell.halfSize / 2.0;
        for (double dx : {-h, h}) {
            for (double dy : {-h, h}) {
                Vec2 c = {cell.center.x + dx, cell.center.y + dy};
                double d = signedDistToPolygon(c, localVertices);
                if (d > bestDist) bestDist = d;
                pq.push({c, h, d});
            }
        }
    }

    // Filter to 8-16 highest
    std::sort(poles.begin(), poles.end(), [](const Pole& a, const Pole& b) { return a.radius > b.radius; });
    if (poles.size() > 16) poles.resize(16);
}

// --- CompositeShape Implementation ---

void CompositeShape::warmUp() const {
    for(const auto& part : parts) part.warmUp();
    // updateCache của CompositeShape (tính totalAABB) đã nằm trong setTransform, 
    // nhưng để chắc chắn ta có thể gọi logic tính AABB ở đây nếu cần.
}

void CompositeShape::addPart(const std::vector<Vec2>& verts) {
    ConvexPolygon p;
    p.setVertices(verts);
    parts.push_back(p);
}

void CompositeShape::setTransform(Vec2 pos, double angle) {
    this->pos = pos;
    this->angle = angle;

    double minX = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();

    for (auto& part : parts) {
        part.setTransform(pos, angle);
        
        const AABB& aabb = part.getAABB();
        if (aabb.min.x < minX) minX = aabb.min.x;
        if (aabb.max.x > maxX) maxX = aabb.max.x;
        if (aabb.min.y < minY) minY = aabb.min.y;
        if (aabb.max.y > maxY) maxY = aabb.max.y;
    }
    totalAABB = {{minX, minY}, {maxX, maxY}};
}

void CompositeShape::precomputeAllPoles(double precision_rel) {
    allPoles.clear();
    diameter = 0.0;
    convexHullArea = 0.0;
    for (auto& part : parts) {
        part.precomputePolesAndMetrics(precision_rel);
        allPoles.insert(allPoles.end(), part.poles.begin(), part.poles.end());
        diameter = std::max(diameter, part.diameter);
        convexHullArea += part.convexHullArea;
    }
}