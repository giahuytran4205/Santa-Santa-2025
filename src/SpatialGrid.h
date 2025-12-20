// src/SpatialGrid.h
#ifndef SPATIAL_GRID_H
#define SPATIAL_GRID_H

#include <vector>
#include <cmath>
#include "Geometry.h"

struct SpatialGrid {
    double cellSize;
    int rows, cols;
    double halfSize;
    // Mảng 1 chiều lưu danh sách các index vật thể trong từng ô
    // grid[idx] -> vector<int> chứa index của items
    std::vector<std::vector<int>> cells;

    SpatialGrid(double containerSize, double maxItemDiameter) {
        // Kích thước mỗi ô lưới ít nhất phải bằng đường kính vật lớn nhất
        // để đảm bảo vật không bao giờ nằm quá 2 ô theo chiều ngang/dọc
        cellSize = maxItemDiameter;
        halfSize = containerSize / 2.0;
        
        // Tính số hàng/cột (thêm padding để an toàn)
        int dim = static_cast<int>(std::ceil(containerSize / cellSize)) + 2;
        rows = dim;
        cols = dim;
        
        cells.resize(rows * cols);
    }

    // Chuyển tọa độ thực sang index ô lưới
    int getCellIndex(Vec2 pos) const {
        // Dời tọa độ từ [-half, half] về [0, size]
        int c = static_cast<int>((pos.x + halfSize) / cellSize) + 1;
        int r = static_cast<int>((pos.y + halfSize) / cellSize) + 1;
        
        // Clamp (kẹp) giá trị để không bị out of bound
        if (c < 0) c = 0; if (c >= cols) c = cols - 1;
        if (r < 0) r = 0; if (r >= rows) r = rows - 1;
        
        return r * cols + c;
    }

    // Thêm vật thể vào lưới
    void insert(int itemIdx, Vec2 pos) {
        int idx = getCellIndex(pos);
        cells[idx].push_back(itemIdx);
    }

    // Lấy danh sách các ứng viên va chạm (candidates)
    // Bao gồm ô chứa pos và 8 ô xung quanh
    void getCandidates(Vec2 pos, std::vector<int>& outCandidates) const {
        int c = static_cast<int>((pos.x + halfSize) / cellSize) + 1;
        int r = static_cast<int>((pos.y + halfSize) / cellSize) + 1;

        // Duyệt 3x3 ô xung quanh
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                int nr = r + dr;
                int nc = c + dc;

                if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                    int idx = nr * cols + nc;
                    const auto& bucket = cells[idx];
                    outCandidates.insert(outCandidates.end(), bucket.begin(), bucket.end());
                }
            }
        }
    }
};

#endif // SPATIAL_GRID_H