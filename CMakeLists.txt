cmake_minimum_required(VERSION 3.3)
project(learnedspatial)

include_directories(include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I. -std=c++17")


add_executable(learnedspatial src/main.cpp)
