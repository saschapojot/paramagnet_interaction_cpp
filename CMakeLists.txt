cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories("/home/cywanag/liuxi/local/include/eigen-3.4.0")
include_directories("/home/cywanag/liuxi/gccNew/include")
include_directories("/home/cywanag/liuxi/msgpack-c-cpp_master/include")
set_source_files_properties(
        dbExchangeModel.cpp
        PROPERTIES
        COMPILE_OPTIONS "-Ofast"
)

add_executable(mc main.cpp dbExchangeModel.cpp dbExchangeModel.hpp)

