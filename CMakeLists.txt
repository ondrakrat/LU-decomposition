cmake_minimum_required(VERSION 3.8)
project(hw1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11")

set(SOURCE_FILES LU_skeleton.cpp)
add_executable(hw1 ${SOURCE_FILES})