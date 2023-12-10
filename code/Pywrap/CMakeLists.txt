
file(GLOB pywrap_cpp ${CMAKE_CURRENT_LIST_DIR}/*.cpp)


target_sources(tmos 
    PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/wrap_base.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/wrap_surface.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/wrap_pili.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/wrap_mdynamics.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/tmos.cpp"
    )

if(${USE_VTK}) 
    target_sources(tmos 
        PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/wrap_vtkwriter.cpp"
        )
else()
    target_sources(tmos 
        PUBLIC
        "${CMAKE_CURRENT_LIST_DIR}/no_vtk.cpp"
        )
endif()

target_include_directories(tmos PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
    
    
