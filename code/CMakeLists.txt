cmake_minimum_required(VERSION 3.5)

project(Pili
    VERSION 0.1.0
    LANGUAGES CXX
)


# find_program(CMAKE_CXX_COMPILER NAMES "g++" "gcc" "x86_64-conda_cos6-linux-gnu-c++" "gcc-8")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS ON) # TODO uses .a, switch to .so
find_package(Boost 1.67 REQUIRED COMPONENTS unit_test_framework)
find_package(NLopt CONFIG REQUIRED)

if(USE_VTK)
    # set(VTK_DIR "${CMAKE_CURRENT_LIST_DIR}/lib/VTK-build" CACHE PATH "VTK CMake package directory")
    find_package(VTK REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type")

# set(USE_VTK false CACHE BOOL "compile with vtk")
#get_filename_component(PILI_DIR ../ ABSOLUTE CACHE)
# get_filename_component(SRC_DIR ${CMAKE_CURRENT_LIST_DIR} ABSOLUTE CACHE)

# set(BOOST_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../lib/boost_1_67_0" CACHE PATH "Boost installation root")
# set(Boost_USE_STATIC_LIBS ON)

# set(NLopt_DIR "${CMAKE_CURRENT_LIST_DIR}/lib/nlopt-2.6.1" CACHE PATH "NLopt CMake package directory")


# add_compile_options(-fPIC)

include_directories(${CMAKE_CURRENT_LIST_DIR}/lib) # eigen
include_directories(${Boost_INCLUDE_DIRS})

set(PYBIND11_CPP_STANDARD -std=c++14)
add_subdirectory(pybind11)

# compile flags
# TODO add optimisation flags?
Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
message(${CMAKE_CXX_FLAGS})

# including a source file here bypasses a "no sources included" error 
pybind11_add_module(tmos SHARED Base/vector3d.cpp)

#https://crascit.com/2016/01/31/enhanced-source-file-handling-with-target_sources/
# include source files for target tmos from all subdirectories
include(Base/CMakeLists.txt)
include(Surface/CMakeLists.txt)
include(Bacterium/CMakeLists.txt)
include(Pywrap/CMakeLists.txt)

if(USE_VTK)
    include(Vtk/CMakeLists.txt)
    target_link_libraries(tmos PRIVATE ${VTK_LIBRARIES})
endif()

# needs to link to external library
target_link_libraries(tmos PRIVATE NLopt::nlopt)

# set_target_properties(tmos PROPERTIES SUFFIX ".so")
set_target_properties(tmos PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Hierarchy of modules
# ?
set(base_inc ${CMAKE_CURRENT_SOURCE_DIR}/Base)
set(surface_inc ${base_inc} ${CMAKE_CURRENT_SOURCE_DIR}/Surface)
set(bacterium_inc ${surface_inc} ${CMAKE_CURRENT_SOURCE_DIR}/Bacterium)
set(vtk_inc ${bacterium_inc} ${CMAKE_CURRENT_SOURCE_DIR}/Vtk)

# TODO
# add_subdirectory(Base)
# add_subdirectory(Surface)
# add_subdirectory(Bacterium)
# add_subdirectory(Pywrap)

# testing 
enable_testing()
add_subdirectory(UnitTest)

# debug
# include(check_vars.cmake)

