cmake_minimum_required(VERSION 3.26)
project(paramagnet_interaction_cpp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)






#SET( CMAKE_CXX_FLAGS "-ggdb")
include_directories("/home/polya/cppLocal/eigen-3.4.0")

include_directories("/home/polya/cppLocal/msgpack-c-cpp_master/include")

#include_directories("/home/polya/cppLocal/capnproto/include")

set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/polya/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

#find_package(Boost  REQUIRED COMPONENTS ALL)
find_package(Boost  REQUIRED serialization filesystem)
if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()


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

set_source_files_properties(
        s2Eigen.cpp
        PROPERTIES
        COMPILE_OPTIONS "-Ofast"
)
add_executable(mc main.cpp dbExchangeModel.cpp)



add_executable(reader reader.cpp sAvg.cpp)


add_executable(band bandUnfoldingAndC.cpp s2Eigen.cpp dbExchangeModel.cpp)

file(GLOB LIBS_A "${Boost_LIBRARY_DIR}/*serialization.a")
file(GLOB LIBS_SO "${Boost_LIBRARY_DIR}/*serialization.so")

target_link_libraries(mc ${LIBS_A})
target_link_libraries(mc ${LIBS_SO})
target_link_libraries(band ${LIBS_A})
target_link_libraries(band ${LIBS_SO})


file(GLOB LIBS_A "${Boost_LIBRARY_DIR}/*filesystem.a")
file(GLOB LIBS_SO "${Boost_LIBRARY_DIR}/*filesystem.so")

target_link_libraries(mc ${LIBS_A})
target_link_libraries(mc ${LIBS_SO})
target_link_libraries(band ${LIBS_A})
target_link_libraries(band ${LIBS_SO})

#
#file(GLOB LIBS_A "/home/polya/cppLocal/capnproto/lib/*.a")
#file(GLOB LIBS_SO "/home/polya/cppLocal/capnproto/lib/*.so")
#
#target_link_libraries(reader ${LIBS_A})
#target_link_libraries(mc ${LIBS_SO})
