

file(GLOB surface_cpp ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
file(GLOB surface_hpp ${CMAKE_CURRENT_LIST_DIR}/*.hpp)


#tmp
#list(REMOVE_ITEM surface_cpp hexsurface.cpp)
#list(REMOVE_ITEM surface_hpp hexsurface.hpp)

target_sources(tmos 
    PRIVATE
    ${surface_cpp}
    PUBLIC
    ${surface_hpp}
    )

target_include_directories(tmos PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
    
    
