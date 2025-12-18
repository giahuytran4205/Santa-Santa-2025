// main.cpp
#include <iostream>
#include "Geometry.h"
#include "Collision.h"
#include "Tree.h"
#include "Sparrow.h"
#include "SVGExport.h"

int main() {
    // 1. Tạo dữ liệu (15 cây thông)
    std::vector<CompositeShape> items;
    for(int i=0; i<5; ++i) items.push_back(createTree());

    // 2. Cấu hình
    SparrowConfig config;
    config.container_size = 1.2; // Kích thước khởi tạo (nên đủ lớn)
    config.max_iter = 1000;       // Iteration cho hàm separate (Inner loop)
    config.n_samples = 2048;
    config.n_threads = 8;
    config.Nc = 20;         // Reduce inner iters
    config.Nx = 10;
    config.Kc = 3;
    config.Kx = 2;
    config.TLx = 10.0;      // 10 seconds for explore
    config.TLc = 10.0;      // 10 seconds for compress
    config.Rx = 0.0005;      // Slightly larger shrink to converge faster
    config.Rs_c = 0.005;
    config.Re_c = 0.0005;

    // 3. Chạy Solver
    SparrowSolver solver(items, config);
    solver.solve(); // Chạy toàn bộ quy trình Algo 10->13
    exportToSVG("../output/ket_qua_sparrow.svg", solver, 40.0);

    // 4. Kết quả nằm trong solver.items và solver.getContainerSize()
    return 0;
}