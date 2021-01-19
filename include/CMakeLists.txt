cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(jacobi_pd)

include_directories()
link_directories()
add_definitions()


add_executable (jacobi_pd main.cpp jacobi_pd.hpp matrix_alloc_jpd.hpp)
target_link_libraries (jacobi_pd)




