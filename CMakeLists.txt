cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home/polya/boost_1_82_0")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home/cywanag/liuxi/local/boost_1_84_0/libs")
FIND_PACKAGE(Boost)

IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

include_directories("/usr/local/include/eigen-3.4.0")
include_directories("/home/polya/cppLocal/msgpack-c-cpp_master/include")
set_source_files_properties(
        dbExchangeModel.cpp
        PROPERTIES
        COMPILE_OPTIONS "-Ofast"
)

add_executable(mc main.cpp dbExchangeModel.cpp dbExchangeModel.hpp)
add_executable(reader reader.cpp sAvg.hpp sAvg.cpp)


