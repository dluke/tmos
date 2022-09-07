
file(GLOB vtk_cpp ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
file(GLOB vtk_hpp ${CMAKE_CURRENT_LIST_DIR}/*.hpp)

target_sources(tmos 
    PRIVATE
    "${CMAKE_CURRENT_LIST_DIR}/vtkwriter.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/vtksurface.cpp"
    PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/vtkwriter.hpp"
    "${CMAKE_CURRENT_LIST_DIR}/vtksurface.hpp"
    )

target_include_directories(tmos PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
    
    
