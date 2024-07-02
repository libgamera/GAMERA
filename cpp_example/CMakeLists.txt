cmake_minimum_required(VERSION 3.17...3.25 FATAL_ERROR)

project(EXAMPLE VERSION 0.1.0 LANGUAGES CXX)

find_package(GAMERA CONFIG REQUIRED)

add_executable(example example.cpp)
target_link_libraries(example PRIVATE gamera::gamera)
