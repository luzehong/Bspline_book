cmake_minimum_required(VERSION 3.22)
project(Bspline_book)

set(CMAKE_CXX_STANDARD 14)
link_directories("/usr/local/include")
add_executable(Bspline_book main.cpp include/class_lis/bsplinesurface.h)