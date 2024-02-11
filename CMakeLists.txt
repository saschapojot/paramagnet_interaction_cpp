cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 20)
include_directories("/home/cywanag/liuxi/local/include/eigen-3.4.0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_executable(mc main.cpp dbExchangeModel.cpp dbExchangeModel.hpp)
#set(BOOST_ROOT /home/cywanag/liuxi/local/boost_1_84_0)

#target_include_directories(mc PUBLIC /home/cywanag/liuxi/local/boost_1_84_0/buildedLibs/include)
#find_package(Boost REQUIRED COMPONENTS system serialization)

#if(Boost_FOUND)
 #   include_directories(${Boost_INCLUDE_DIRS})
  #  add_executable(mc main.cpp dbExchangeModel.cpp dbExchangeModel.hpp)
   # target_link_directories(mc PUBLIC /home/cywanag/liuxi/local/boost_1_84_0/buildedLibs/lib)
#else()
 #   message(FATAL_ERROR "Boost not found. Please check your Boost installation.")
#endif()