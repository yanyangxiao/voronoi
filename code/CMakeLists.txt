cmake_minimum_required (VERSION 3.1)

project("voronoi")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-D_USE_MATH_DEFINES)

#BOOST
find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDES})

#CGAL
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDES})
link_directories(${CGAL_LIBRARIES})
link_libraries(${CGAL_LIBRARIES})

find_package(OpenMP REQUIRED)

add_executable(voronoi2-demo voronoi2-demo.cpp)
target_link_libraries(voronoi2-demo  ${CGAL_LIBRARIES})
target_link_libraries(voronoi2-demo OpenMP::OpenMP_CXX)

add_executable(voronoi3-demo voronoi3-demo.cpp)
target_link_libraries(voronoi3-demo  ${CGAL_LIBRARIES})
target_link_libraries(voronoi3-demo OpenMP::OpenMP_CXX)

add_executable(power2-demo power2-demo.cpp)
target_link_libraries(power2-demo  ${CGAL_LIBRARIES})
target_link_libraries(power2-demo OpenMP::OpenMP_CXX)

add_executable(power3-demo power3-demo.cpp)
target_link_libraries(power3-demo  ${CGAL_LIBRARIES})
target_link_libraries(power3-demo OpenMP::OpenMP_CXX)

### OPENMESH
include_directories("D:/library/OpenMesh/OpenMesh9/src/")
link_directories("D:/library/OpenMesh/OpenMesh9/lib/")
add_definitions(-D_USE_MATH_DEFINES)### OPENMESH
include_directories("D:/library/OpenMesh/OpenMesh9/src/")
link_directories("D:/library/OpenMesh/OpenMesh9/lib/")
add_definitions(-D_USE_MATH_DEFINES)

add_executable(voronois-demo voronois-demo.cpp)
target_link_libraries(voronois-demo  ${CGAL_LIBRARIES})
target_link_libraries(voronois-demo debug OpenMeshCored.lib optimized OpenMeshCore.lib)
target_link_libraries(voronois-demo OpenMP::OpenMP_CXX)

add_executable(powers-demo powers-demo.cpp)
target_link_libraries(powers-demo  ${CGAL_LIBRARIES})
target_link_libraries(powers-demo debug OpenMeshCored.lib optimized OpenMeshCore.lib)
target_link_libraries(powers-demo OpenMP::OpenMP_CXX)

