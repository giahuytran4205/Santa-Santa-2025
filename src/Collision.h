// Collision.h
#ifndef COLLISION_H
#define COLLISION_H

#include "Geometry.h"

// Kiểm tra SAT giữa 2 đa giác lồi
bool checkSAT(const ConvexPolygon& polyA, const ConvexPolygon& polyB);

// Kiểm tra va chạm giữa 2 vật thể phức hợp
bool checkCollisionComposite(const CompositeShape& shapeA, const CompositeShape& shapeB);

// THÊM/THAY THẾ: Hybrid quantification
double quantify_collision(const CompositeShape& shapeA, const CompositeShape& shapeB);

#endif // COLLISION_H