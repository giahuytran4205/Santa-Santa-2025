#ifndef COLLISION_H
#define COLLISION_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "Geometry.h"

// --- Cấu trúc SoA (Structure of Arrays) ---
// Định nghĩa lại CirclesSoA ở đây hoặc dùng forward declaration nếu đã có trong Geometry.h.
// Tuy nhiên, để file này độc lập và dễ include, ta định nghĩa đầy đủ (giống như trong Geometry.h).
// Lưu ý: Nếu compiler báo lỗi redefinition, hãy dùng #ifndef hoặc tách struct này ra file riêng (VD: Common.h).
// Trong ngữ cảnh này, vì Geometry.h đã được include, ta sẽ sử dụng định nghĩa từ Geometry.h
// nếu Geometry.h đã có. Nhưng để an toàn và nhất quán với các file trước, 
// ta giả định Geometry.h CHƯA định nghĩa struct này ở global scope mà chỉ dùng trong CompositeShape.
//
// TỐT NHẤT: Ta sẽ tái sử dụng định nghĩa `CirclesSoA` đã có trong `Geometry.h` 
// (đã được bạn yêu cầu sửa ở bước trước).
// Vì vậy ở đây ta KHÔNG định nghĩa lại `CirclesSoA` nữa mà chỉ sử dụng nó.

// Hàm tính va chạm tối ưu (Fast Circle Physics)
// Đầu vào là hai tập hợp hình tròn (SoA) và ngưỡng epsilon để làm mềm va chạm
// Hàm này được cài đặt trong Collision.cpp
double evaluate_circles_fast(const CirclesSoA& setA, const CirclesSoA& setB, double epsilon);

// Hàm kiểm tra va chạm đa giác truyền thống
// Giữ lại để đảm bảo tính tương thích hoặc dùng cho giai đoạn kiểm tra cuối cùng (nếu cần độ chính xác tuyệt đối)
bool checkCollisionComposite(const CompositeShape& shapeA, const CompositeShape& shapeB);

// Hàm tính toán chi tiết độ chồng lấn giữa 2 shape (dùng logic cũ hoặc wrap logic mới)
double quantify_collision(const CompositeShape& shapeA, const CompositeShape& shapeB);

#endif // COLLISION_H