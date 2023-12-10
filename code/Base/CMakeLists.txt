

file(GLOB base_cpp ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
file(GLOB base_hpp ${CMAKE_CURRENT_LIST_DIR}/*.hpp)

target_sources(tmos 
    PRIVATE
    ${base_cpp}
    PUBLIC
    ${base_hpp}
    )

target_include_directories(tmos PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
    
    
