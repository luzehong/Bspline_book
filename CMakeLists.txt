cmake_minimum_required(VERSION 3.20)
project(Bspline_book)

set(CMAKE_CXX_STANDARD 14)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
add_executable(Bspline_book main.cpp include/class_lis/bsplinesurface.h)