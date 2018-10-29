cmake_minimum_required(VERSION 3.12)
project(High_Performance_Homework)

set(CMAKE_CXX_STANDARD 14)
include_directories(/share/apps/openmpi/include)
SET(CMAKE_CXX_COMPILER mpicxx)
add_executable(mat_pro exp_2/code/A1.cpp)