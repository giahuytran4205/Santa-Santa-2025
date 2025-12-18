// SVGExport.cpp
#include "SVGExport.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

// Helper: Tính diện tích đa giác (Shoelace formula)
// double getPolygonArea(const std::vector<Vec2>& points) {
//     double area = 0.0;
//     for (size_t i = 0; i < points.size(); ++i) {
//         size_t j = (i + 1) % points.size();
//         area += points[i].x * points[j].y;
//         area -= points[j].x * points[i].y;
//     }
//     return std::abs(area) / 2.0;
// }

// Helper: Format số thực (ví dụ 3 chữ số thập phân)
std::string fmt(double val, int precision = 3) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << val;
    return ss.str();
}

void exportToSVG(const std::string& filename, const SparrowSolver& solver, double scale) {
    const auto& items = solver.getItems();
    double containerSize = solver.getContainerSize();
    
    // 1. Tính toán thống kê (Header Info)
    double totalItemArea = 0.0;
    for (const auto& item : items) {
        for (const auto& part : item.parts) {
            totalItemArea += getPolygonArea(part.getLocalVertices());
        }
    }
    double containerArea = containerSize * containerSize;
    double density = (totalItemArea / containerArea) * 100.0;

    // 2. Thiết lập ViewBox
    // Container của Solver là [-size/2, size/2]. 
    // Ta thêm margin phía trên để chứa text header.
    double halfSize = containerSize / 2.0;
    double margin = containerSize * 0.1; // 10% lề
    double headerSpace = containerSize * 0.15; // Chỗ cho text

    double viewBoxMinX = -halfSize - margin;
    double viewBoxMinY = -halfSize - margin - headerSpace; // Mở rộng lên trên
    double viewBoxWidth = containerSize + margin * 2;
    double viewBoxHeight = containerSize + margin * 2 + headerSpace;

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening " << filename << std::endl;
        return;
    }

    // --- Bắt đầu ghi SVG ---
    file << "<svg viewBox='" << fmt(viewBoxMinX) << " " << fmt(viewBoxMinY) << " " 
         << fmt(viewBoxWidth) << " " << fmt(viewBoxHeight) << "' "
         << "xmlns='http://www.w3.org/2000/svg'>" << std::endl;

    // A. Header Info Text
    file << "<text font-family='monospace' font-size='" << fmt(containerSize * 0.04) << "' "
         << "font-weight='500' x='" << fmt(-halfSize) << "' y='" << fmt(-halfSize - containerSize * 0.05) << "'>"
         << "h: " << fmt(containerSize) << " | w: " << fmt(containerSize) 
         << " | d: " << fmt(density) << "% | sparrow_result"
         << "</text>" << std::endl;

    // B. Container Group
    file << "<g id='container_0'>" << std::endl;
    file << "  <rect x='" << fmt(-halfSize) << "' y='" << fmt(-halfSize) << "' "
         << "width='" << fmt(containerSize) << "' height='" << fmt(containerSize) << "' "
         << "fill='#D3D3D3' stroke='black' stroke-width='" << fmt(containerSize * 0.005) << "'/>" << std::endl;
    file << "  <title>container, size: " << fmt(containerSize) << "</title>" << std::endl;
    file << "</g>" << std::endl;

    // C. Items Group
    file << "<g id='items'>" << std::endl;
    
    // C.1 Definitions (<defs>)
    // Định nghĩa hình dáng gốc của từng item (ở local space)
    file << "  <defs>" << std::endl;
    for (size_t i = 0; i < items.size(); ++i) {
        file << "    <g id='def_item_" << i << "'>" << std::endl;
        // Gom tất cả các phần (parts) của item thành 1 path hoặc group
        for (const auto& part : items[i].parts) {
            file << "      <path d='M";
            const auto& verts = part.getLocalVertices();
            for (size_t k = 0; k < verts.size(); ++k) {
                file << fmt(verts[k].x) << "," << fmt(verts[k].y);
                if (k < verts.size() - 1) file << " L";
            }
            file << " z' fill='#7A7A7A' fill-opacity='0.5' fill-rule='nonzero' "
                 << "stroke='black' stroke-width='" << fmt(containerSize * 0.002) << "'/>" << std::endl;
        }
        file << "    </g>" << std::endl;
    }
    file << "  </defs>" << std::endl;

    // C.2 Instantiation (<use>)
    // Đặt item vào vị trí thực tế bằng transform
    for (size_t i = 0; i < items.size(); ++i) {
        const auto& item = items[i];
        
        // Chuyển đổi góc từ radian sang độ
        double angleDeg = item.angle * (180.0 / PI);
        
        file << "  <use href='#def_item_" << i << "' "
             << "transform='translate(" << fmt(item.pos.x) << " " << fmt(item.pos.y) << "), "
             << "rotate(" << fmt(angleDeg) << ")'>" << std::endl;
        
        file << "    <title>item, id: " << i << ", transf: [r: " << fmt(angleDeg) 
             << "°, t: (" << fmt(item.pos.x) << ", " << fmt(item.pos.y) << ")]</title>" << std::endl;
        file << "  </use>" << std::endl;
    }
    file << "</g>" << std::endl;

    // Kết thúc file
    file << "</svg>" << std::endl;
    file.close();
    
    std::cout << "Exported SVG: " << filename << std::endl;
}