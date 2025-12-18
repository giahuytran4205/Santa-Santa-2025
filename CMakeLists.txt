cmake_minimum_required(VERSION 3.12)
project(SparrowProject)

# 1. Cấu hình C++ Standard (Dùng C++17 hoặc 20)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2. Tìm và kích hoạt OpenMP
find_package(OpenMP REQUIRED)

# 3. TỰ ĐỘNG QUÉT FILE SOURCE (Đây là phần bạn cần)
# GLOB_RECURSE: Quét cả các thư mục con
# CONFIGURE_DEPENDS: Tự động chạy lại cmake khi bạn thêm/xóa file (CMake 3.12+)
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.h")

# 4. Tạo file thực thi
add_executable(sparrow_run ${SOURCES} ${HEADERS})

# 5. Link thư viện OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(sparrow_run PUBLIC OpenMP::OpenMP_CXX)
endif()

# Tối ưu hóa (tương đương -O3)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()