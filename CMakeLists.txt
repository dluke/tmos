cmake_minimum_required(VERSION 3.2.2)

find_program(CMAKE_CXX_COMPILER NAMES "g++" "gcc" "x86_64-conda_cos6-linux-gnu-c++" "gcc-8")
set(CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_PREFIX "/homes/dbarton/miniconda3" CACHE PATH "?")

project(Pili)

find_package (Python3 3.8 COMPONENTS Interpreter Development)
# set(USE_PYTHON_INCLUDE_DIR ${Python3_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type")

set(USE_VTK false CACHE BOOL "compile with vtk")
#get_filename_component(PILI_DIR ../ ABSOLUTE CACHE)
get_filename_component(SRC_DIR ${CMAKE_CURRENT_LIST_DIR} ABSOLUTE CACHE)

find_package(Boost)
find_library(Boost_TEST_LIB
    NAMES boost_unit_test_framework
    PATHS /usr/lib/x86_64-linux-gnu/
)

if(${USE_VTK}) 
    SET(VTK_DIR ${CMAKE_CURRENT_LIST_DIR}/../lib/VTK-build) 
    find_package(VTK REQUIRED PATHS ${CMAKE_CURRENT_LIST_DIR}/../lib/VTK-build)
endif()

add_compile_options(-fPIC)

include_directories(/usr/include/python3.8)
include_directories($ENV{CONDA_PREFIX}/include/python3.8)
include_directories(${CMAKE_CURRENT_LIST_DIR}/lib) # eigen
include_directories(${Boost_INCLUDE_DIR})
# include_directories(${EIGEN3_INCLUDE_DIR})

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

if(${USE_VTK})
    include(Vtk/CMakeLists.txt)
    target_link_libraries(tmos PRIVATE ${VTK_LIBRARIES})
endif()

# needs to link to external library
find_library(NLOPTLIB
    NAMES nlopt
    PATHS /usr/local/lib /usr/lib)
message(${NLOPTLIB})
target_link_libraries(tmos PRIVATE ${NLOPTLIB})

set_target_properties(tmos PROPERTIES SUFFIX ".so")

# Hierarchy of modules
# ?
set(base_inc ${CMAKE_CURRENT_SOURCE_DIR}/Base)
set(surface_inc ${base_inc} ${CMAKE_CURRENT_SOURCE_DIR}/Surface)
set(bacterium_inc ${surface_inc} ${CMAKE_CURRENT_SOURCE_DIR}/Bacterium)
set(vtk_inc ${bacterium_inc} ${CMAKE_CURRENT_SOURCE_DIR}/Vtk)

# testing 
# enable_testing()

add_subdirectory(${PROJECT_SOURCE_DIR}/UnitTest)


#add_custom_target(build_tests)
#add_dependencies(build_tests ${PROJECT_SOURCE_DIR}/UnitTest)

# debug
# include(check_vars.cmake)

