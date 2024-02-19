cmake_minimum_required(VERSION 3.26)
project(pilocal)

set(CMAKE_CXX_STANDARD 20)
include_directories("/home/polya/cppLocal/msgpack-c-cpp_master/include")
include_directories("/usr/local/include/eigen-3.4.0")


set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIRS "${Boost_ROOT}/include")
    set(Boost_LIBRARY_DIRS "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home/cywanag/liuxi/local/boost_1_84_0/libs")


find_package(Boost  REQUIRED COMPONENTS ALL)
if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIRS}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIRS}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
else()
    message(WARNING "boost not found.")
endif()


#link_directories("/home/polya/cppLocal/boostBuild/lib")
#file(GLOB BoostLib_a "/home/polya/cppLocal/boostBuild/lib/*.a")
#file(GLOB BoostLib_so "/home/polya/cppLocal/boostBuild/lib/*.so")
add_executable(mc main.cpp
        MyClass.hpp
        MyClass.cpp
#        ${BoostLib_a}
#        ${BoostLib_so}
)
#add_executable(mc main.cpp)
file(GLOB LIBS_TO_LINK_SO "/home/polya/cppLocal/boostBuild/lib/*.so")
foreach(LIB ${LIBS_TO_LINK_SO})
    target_link_libraries(mc PUBLIC ${LIB})
endforeach()
#
file(GLOB LIBS_TO_LINK_A "/home/polya/cppLocal/boostBuild/lib/*.a")
foreach(LIB ${LIBS_TO_LINK_A})
    target_link_libraries(mc PUBLIC ${LIB})
endforeach()
#


#add_executable(read read.cpp)
