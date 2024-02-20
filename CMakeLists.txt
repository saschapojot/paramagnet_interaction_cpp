cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_STANDARD_REQUIRED ON)


#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home/polya/cppLocal/boost_1_82_0/buildedLibs/include")
#
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home/polya/cppLocal/boost_1_82_0/buildedLibs/lib")
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home/polya/cppLocal/boost_1_82_0")


FIND_PACKAGE(Boost)

IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

#SET( CMAKE_CXX_FLAGS "-ggdb")
include_directories("/home/polya/cppLocal/eigen-3.4.0")
set(GSL_INCLUDE_DIR "/home/polya/cppLocal/GSL/include")
set(GSL_LIBRARY_DIR "/home/polya/cppLocal/GSL/lib")



include_directories("/home/polya/cppLocal/msgpack-c-cpp_master/include")
#include_directories("/home/polya/cppLocal/protobuf/include")



set_source_files_properties(
        dbExchangeModel.cpp
        PROPERTIES
        COMPILE_OPTIONS "-Ofast"
)

set_source_files_properties(
        sAvg.cpp
        PROPERTIES
        COMPILE_OPTIONS "-Ofast"
)

link_directories("/home/polya/cppLocal/protobuf/lib")
add_executable(mc main.cpp dbExchangeModel.cpp dbExchangeModel.hpp)


add_executable(reader reader.cpp sAvg.cpp sAvg.hpp dbExchangeModel.hpp)

#add_executable(tut protoTut.cpp addressbook.pb.h addressbook.pb.cc)


target_include_directories(mc PRIVATE ${GSL_INCLUDE_DIR})

# Link against GSL libraries
target_link_directories(mc PRIVATE ${GSL_LIBRARY_DIR})
target_link_libraries(mc PRIVATE gsl gslcblas m)


#file(GLOB LIB_FILES_A "/home/polya/cppLocal/protobuf/lib/*.a")
#
#foreach(LIB_FILE ${LIB_FILES_A})
#    target_link_libraries(tut ${LIB_FILE})
#endforeach()
#
#file(GLOB LIB_FILES_SO "/home/polya/cppLocal/protobuf/lib/*.so")
#
#foreach(LIB_FILE ${LIB_FILES_SO})
#    target_link_libraries(tut ${LIB_FILE})
#endforeach()

#file(GLOB LIB_FILES_LA "/home/polya/cppLocal/protobuf/lib/*.la")
#
#foreach(LIB_FILE ${LIB_FILES_LA})
#    target_link_libraries(tut ${LIB_FILE})
#endforeach()