// MathUtils.h
#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>
#include <algorithm>

constexpr double PI = M_PI;

// Macro ép buộc nội tuyến (Force Inline)
#if defined(_MSC_VER)
    #define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define FORCE_INLINE inline __attribute__((always_inline))
#else
    #define FORCE_INLINE inline
#endif

// --- Vec2 ---
struct Vec2 {
    double x, y;

    FORCE_INLINE Vec2 operator+(const Vec2& v) const { return {x + v.x, y + v.y}; }
    FORCE_INLINE Vec2 operator-(const Vec2& v) const { return {x - v.x, y - v.y}; }
    FORCE_INLINE Vec2 operator*(double s) const { return {x * s, y * s}; }
    FORCE_INLINE double dot(const Vec2& v) const { return x * v.x + y * v.y; }
    
    FORCE_INLINE Vec2 rotate(double sin_a, double cos_a) const {
        return {x * cos_a - y * sin_a, x * sin_a + y * cos_a};
    }

    FORCE_INLINE double length() const { return std::sqrt(x * x + y * y); }  // THÊM

    bool operator==(const Vec2& other) const {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Vec2& other) const {
        return x != other.x || y != other.y;
    }

    bool operator<(const Vec2& other) const {
        if (x != other.x) return x < other.x;
        return y < other.y;
    }
};

// --- AABB ---
struct AABB {
    Vec2 min, max;
    
    FORCE_INLINE bool overlaps(const AABB& other) const {
        return (min.x <= other.max.x && max.x >= other.min.x &&
                min.y <= other.max.y && max.y >= other.min.y);
    }
};

// --- Transform ---
struct Transform {
    Vec2 pos;
    double angle;
    double sin_a, cos_a;

    Transform(Vec2 p = {0.0, 0.0}, double a = 0.0) : pos(p), angle(a) {
        updateTrig();
    }

    void set(Vec2 p, double a) {
        if (p.x != pos.x || p.y != pos.y || a != angle) {
            pos = p;
            if (angle != a) {
                angle = a;
                updateTrig();
            }
        }
    }
    
    void updateTrig() {
        sin_a = std::sin(angle);
        cos_a = std::cos(angle);
    }

    bool operator==(const Transform& other) const {
        return pos == other.pos && angle == other.angle;
    }

    bool operator<(const Transform& other) const {
        if (pos != other.pos) return pos < other.pos;
        return angle < other.angle;
    }
};

#endif // MATH_UTILS_H