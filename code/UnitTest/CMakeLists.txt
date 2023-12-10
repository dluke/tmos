

# extra


add_executable(tboost
    test_boost_load.cpp
    )

target_include_directories(tboost PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(tboost PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)

add_test(NAME boost COMMAND tboost)


#0

add_executable(tbase
    test_base.cpp
    ${base_cpp}
    )

target_include_directories(tbase PRIVATE ${base_inc})
target_link_libraries(tbase PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)

add_test(NAME base COMMAND tbase)


add_executable(tsurface
    test_surface.cpp
    ${base_cpp}
    ${surface_cpp}
    ${bacterium_cpp}
    )

target_include_directories(tsurface PUBLIC ${bacterium_inc} ${PYBIND11_INCLUDE_DIR})
#target_include_directories(tsurface PRIVATE ${surface_inc} ${Boost_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIR})
#target_include_directories(tsurface PRIVATE ../Surface/shape.cpp ${Boost_INCLUDE_DIRS})
target_link_libraries(tsurface PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)
target_link_libraries(tsurface PRIVATE ${NLOPTLIB})

add_test(NAME surface COMMAND tsurface)


#1

include_directories(../)
include_directories(${PILI_DIR}/lib)
include_directories(${PYBIND11_INCLUDE_DIR})

add_executable(tgeometry 
    test_geometry.cpp 
    ${base_cpp}
    ${surface_cpp}
    ${bacterium_cpp}
    )

target_include_directories(tgeometry PRIVATE ${bacterium_inc})
target_link_libraries(tgeometry PRIVATE ${NLOPTLIB})

add_test(NAME geometry COMMAND tgeometry)

#

add_executable(tdetach
    test_detach.cpp 
    ${base_cpp}
    ${surface_cpp}
    ${bacterium_cpp}
    )

target_include_directories(tdetach PRIVATE ${bacterium_inc})
target_link_libraries(tdetach PRIVATE ${NLOPTLIB})

add_test(NAME detach COMMAND tdetach)


# 1.1.

add_executable(tcell
    test_cell.cpp
    ${surface_cpp}
    ${base_cpp}
    ${bacterium_cpp}
    )

target_include_directories(tcell PRIVATE ${bacterium_inc} ${PYBIND11_INCLUDE_DIR})
target_link_libraries(tcell PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)

add_test(NAME cell COMMAND tcell)

# 1.2.

add_executable(tevent
    test_event.cpp
    ${base_cpp}
    ${surface_cpp}
    ${bacterium_cpp}
    )

target_include_directories(tevent PRIVATE ${bacterium_inc} ${PYBIND11_INCLUDE_DIR})
target_link_libraries(tevent PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)

add_test(NAME event COMMAND tevent)




#2

add_executable(tstr 
    test_str.cpp 
    ${base_cpp}
    ${surface_cpp}
    ${bacterium_cpp}
    )

target_include_directories(tstr PRIVATE ${bacterium_inc} ${PYBIND11_INCLUDE_DIR})
target_link_libraries(tstr PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)

add_test(NAME stroutput COMMAND tstr)

#3 

add_executable(tfunction
    test_newton.cpp
    )

add_test(NAME functions COMMAND tfunction)

#4


add_executable(tvtk
    test_vtk.cpp
    ${base_cpp}
    ${surface_cpp}
    ${bacterium_cpp}
    ../Vtk/vtksurface.cpp
    )

target_include_directories(tvtk PUBLIC ${vtk_inc} ${Boost_INCLUDE_DIRS})
target_link_libraries(tvtk ${VTK_LIBRARIES} ${Boost_LIBRARIES})

# add_test(NAME vtkmain COMMAND tvtk)


#5
#eigen

pybind11_add_module(teigen SHARED test_eigen.cpp)
set_target_properties(teigen PROPERTIES SUFFIX ".so")
#target_sources(teigen PUBLIC test_eigen.cpp)

#add_test(NAME teigen COMMAND teigen)

#target_include_directories(tvtk PUBLIC ${vtk_inc} ${Boost_INCLUDE_DIRS})
#target_link_libraries(teigen PRIVATE Eigen3::Eigen)




# practice

add_executable(tmap
    test_map.cpp
    )

target_link_libraries(tmap PRIVATE ${Boost_LIBRARY_DIR}/libboost_unit_test_framework.a)

add_test(NAME map COMMAND tmap)



