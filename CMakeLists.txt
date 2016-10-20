cmake_minimum_required(VERSION 2.6)

project(mc CXX)
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -O3")

add_library(mc_ising STATIC mc_ising.cpp)

add_executable(mc.out mc.cpp)
target_link_libraries(mc.out mc_ising)

add_executable(mc_dist.out mc_dist.cpp)
target_link_libraries(mc_dist.out mc_ising)
