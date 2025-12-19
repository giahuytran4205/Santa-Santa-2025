// main.cpp
#include <iostream>
#include "Geometry.h"
#include "Collision.h"
#include "Tree.h"
#include "Sparrow.h"
#include "SVGExport.h"
#include <iomanip>
#include "CLI11.hpp"


int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();
    CLI::App app{"Sparrow Packing Solver"};

    std::string configFilePath = "../config.txt";
    std::string outputFilePath = "../output/solution.svg";
    int num_tree = 3;

    // Định nghĩa các option
    app.add_option("-c,--config", configFilePath, "Path to configuration file");
    app.add_option("-o,--output", outputFilePath, "Path to output SVG file");
    app.add_option("-n,--num-tree", num_tree, "Number of trees");

    // Parse (Macro này tự bắt lỗi và hiện help nếu cần)
    CLI11_PARSE(app, argc, argv);

    // 1. Tạo dữ liệu (15 cây thông)
    std::vector<CompositeShape> items;
    for(int i=0; i < num_tree; ++i) items.push_back(createTree());

    // 2. Cấu hình
    SparrowConfig config;
    config.loadFromFile(configFilePath);

    SparrowSolver solver(items, config);
    solver.solve();
    double side = solver.getContainerSize();
    double score = side * side / items.size();
    std::cout << ">>> Final Score: " << std::setprecision(10) << std::fixed << score << std::endl;
    exportToSVG(outputFilePath, solver, 40.0);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Running time: " << elapsed.count() / 60 << " minutes" << std::endl;

    return 0;
}