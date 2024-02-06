cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 17)

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(BOOST_ROOT /home/polya/boost_1_82_0)
find_package(Boost 1.45.0)
find_package (Threads)
include_directories("/home/polya/eigen-3.4.0")


if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(mc main.cpp
            dbExchangeModel.hpp
            dbExchangeModel.cpp)
    target_link_libraries(mc ${Boost_LIBRARIES})
    target_link_libraries(mc Threads::Threads)
endif()
