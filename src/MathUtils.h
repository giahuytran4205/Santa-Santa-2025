#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>
#include <algorithm>
#include <iostream>

// --- Constants ---
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double INF = 1e30; // Giá trị vô cực đủ lớn

// --- Vec2 Structure ---
struct Vec2 {
    double x, y;

    // Constructors
    Vec2() : x(0), y(0) {}
    Vec2(double _x, double _y) : x(_x), y(_y) {}

    // Operators
    Vec2 operator+(const Vec2& other) const { return {x + other.x, y + other.y}; }
    Vec2 operator-(const Vec2& other) const { return {x - other.x, y - other.y}; }
    Vec2 operator*(double s) const { return {x * s, y * s}; }
    Vec2 operator/(double s) const { return {x / s, y / s}; }
    
    Vec2& operator+=(const Vec2& other) { x += other.x; y += other.y; return *this; }
    Vec2& operator-=(const Vec2& other) { x -= other.x; y -= other.y; return *this; }
    Vec2& operator*=(double s) { x *= s; y *= s; return *this; }

    bool operator==(const Vec2& other) const { return x == other.x && y == other.y; }
    bool operator!=(const Vec2& other) const { return !(*this == other); }
    
    // Lexicographical compare for sorting/map keys
    bool operator<(const Vec2& other) const {
        if (x != other.x) return x < other.x;
        return y < other.y;
    }

    // Methods
    double dot(const Vec2& other) const { return x * other.x + y * other.y; }
    double cross(const Vec2& other) const { return x * other.y - y * other.x; }
    
    double lengthSq() const { return x * x + y * y; }
    double length() const { return std::sqrt(lengthSq()); }
    
    double distSq(const Vec2& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return dx * dx + dy * dy;
    }
    
    double dist(const Vec2& other) const { return std::sqrt(distSq(other)); }

    Vec2 normalized() const {
        double l = length();
        return (l > 1e-9) ? Vec2(x / l, y / l) : Vec2(0, 0);
    }
    
    // Trả về vector pháp tuyến (xoay 90 độ và chuẩn hóa)
    // Dùng để tính normal của cạnh đa giác
    Vec2 normal() const {
        return Vec2(-y, x).normalized();
    }
    
    Vec2 rotated(double angle) const {
        double c = std::cos(angle);
        double s = std::sin(angle);
        return {x * c - y * s, x * s + y * c};
    }
};

// --- AABB (Axis Aligned Bounding Box) ---
struct AABB {
    Vec2 min;
    Vec2 max;

    AABB() : min({INF, INF}), max({-INF, -INF}) {}
    AABB(Vec2 _min, Vec2 _max) : min(_min), max(_max) {}

    double width() const { return max.x - min.x; }
    double height() const { return max.y - min.y; }
    Vec2 center() const { return (min + max) * 0.5; }
    double area() const { return width() * height(); }

    bool overlaps(const AABB& other) const {
        if (max.x < other.min.x || min.x > other.max.x) return false;
        if (max.y < other.min.y || min.y > other.max.y) return false;
        return true;
    }

    bool contains(const Vec2& p) const {
        return p.x >= min.x && p.x <= max.x &&
               p.y >= min.y && p.y <= max.y;
    }
    
    // Mở rộng AABB để chứa thêm điểm p
    void expand(const Vec2& p) {
        if (p.x < min.x) min.x = p.x;
        if (p.y < min.y) min.y = p.y;
        if (p.x > max.x) max.x = p.x;
        if (p.y > max.y) max.y = p.y;
    }
    
    // Mở rộng AABB để chứa thêm AABB khác
    void expand(const AABB& other) {
        if (other.min.x < min.x) min.x = other.min.x;
        if (other.min.y < min.y) min.y = other.min.y;
        if (other.max.x > max.x) max.x = other.max.x;
        if (other.max.y > max.y) max.y = other.max.y;
    }
};

// --- Transform ---
struct Transform {
    Vec2 pos;
    double angle;

    Transform() : pos({0, 0}), angle(0.0) {}
    Transform(Vec2 p, double a) : pos(p), angle(a) {}

    bool operator==(const Transform& other) const {
        return pos == other.pos && std::abs(angle - other.angle) < 1e-9;
    }
};

// --- Helper Functions ---
inline double toRadians(double deg) { return deg * PI / 180.0; }
inline double toDegrees(double rad) { return rad * 180.0 / PI; }

// Utility để in Vec2 ra stream (debug)
inline std::ostream& operator<<(std::ostream& os, const Vec2& v) {
    return os << "(" << v.x << ", " << v.y << ")";
}

#endif // MATH_UTILS_H