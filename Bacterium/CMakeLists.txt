
file(GLOB bacterium_cpp ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
file(GLOB bacterium_hpp ${CMAKE_CURRENT_LIST_DIR}/*.hpp)

target_sources(tmos 
    PRIVATE
    ${bacterium_cpp}
    PUBLIC
    ${bacterium_hpp}
    )

target_include_directories(tmos PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
    
    
