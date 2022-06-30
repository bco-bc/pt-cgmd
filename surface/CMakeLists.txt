cmake_minimum_required(VERSION 3.20)
project(surface)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

include_directories(include)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)

add_library(surface SHARED
            src/vertex.cpp src/edge.cpp src/face.cpp src/triangle.cpp src/polyhedron.cpp
            src/dotted-surface-generator.cpp src/polyhedron-generator.cpp
            src/triangulation.cpp src/nsc.cpp src/pentagon.cpp src/sphere.cpp)

add_subdirectory(tests)

install(TARGETS surface DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
install(DIRECTORY resources DESTINATION ${CMAKE_INSTALL_PREFIX})
