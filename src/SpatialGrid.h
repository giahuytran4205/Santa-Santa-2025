// src/SpatialGrid.h
#ifndef SPATIAL_GRID_H
#define SPATIAL_GRID_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "Geometry.h"

struct SpatialGrid {
    double cellSize;
    int rows, cols;
    double halfSize;
    std::vector<std::vector<int>> cells;

    SpatialGrid(double containerSize, double maxItemDiameter) {
        // [QUAN TRỌNG] Cell size không được nhỏ hơn đường kính vật
        // Nếu maxItemDiameter lỗi (=0), ta fallback về 0.1 để tránh crash
        cellSize = std::max(maxItemDiameter, 0.1); 
        halfSize = containerSize / 2.0;
        
        // [AN TOÀN] Thêm padding rộng (+6 ô) để bao trọn các vật bị đẩy ra ngoài biên
        int dim = static_cast<int>(std::ceil(containerSize / cellSize)) + 6; 
        rows = dim;
        cols = dim;
        
        cells.resize(rows * cols);
    }

    // Chuyển tọa độ thực sang index ô lưới
    int getCellIndex(Vec2 pos) const {
        // Dịch chuyển gốc tọa độ và cộng padding (+3)
        int c = static_cast<int>((pos.x + halfSize) / cellSize) + 3;
        int r = static_cast<int>((pos.y + halfSize) / cellSize) + 3;
        
        // Kẹp giá trị (Clamp) để không bao giờ truy cập ngoài mảng
        if (c < 0) c = 0; else if (c >= cols) c = cols - 1;
        if (r < 0) r = 0; else if (r >= rows) r = rows - 1;
        
        return r * cols + c;
    }

    void insert(int itemIdx, Vec2 pos) {
        cells[getCellIndex(pos)].push_back(itemIdx);
    }

    // Lấy danh sách các vật thể ở ô hiện tại và 8 ô xung quanh
    void getCandidates(Vec2 pos, std::vector<int>& outCandidates) const {
        outCandidates.clear();
        int idx = getCellIndex(pos);
        int c = idx % cols;
        int r = idx / cols;

        // Duyệt 3x3 ô xung quanh
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                int nr = r + dr;
                int nc = c + dc;
                // Kiểm tra biên
                if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
                    const auto& bucket = cells[nr * cols + nc];
                    // Copy dữ liệu
                    outCandidates.insert(outCandidates.end(), bucket.begin(), bucket.end());
                }
            }
        }
    }
};

#endif // SPATIAL_GRID_H