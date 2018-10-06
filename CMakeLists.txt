cmake_minimum_required(VERSION 3.10)
project(making_lattice_test)

set(CMAKE_CXX_STANDARD 11)

add_executable(making_lattice_test main.cpp grcv2d.cpp grcv2d.h naca.cpp naca.h)