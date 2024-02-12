cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories("/home/cywanag/liuxi/local/include/eigen-3.4.0")
include_directories("/home/cywanag/liuxi/gccNew/include")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(mc main.cpp dbExchangeModel.cpp dbExchangeModel.hpp cxx_std_20)

