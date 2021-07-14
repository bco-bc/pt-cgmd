cmake_minimum_required(VERSION 3.19)
project(particles)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /home/andre)

include_directories(${CMAKE_INSTALL_PREFIX}/include)

include_directories(include)
add_library(particles SHARED
        src/atom.cpp src/atomistic.cpp src/bead.cpp src/coarse-grained.cpp src/continuous-protonatable-bead.cpp
        src/discrete-protonatable-bead.cpp src/particle.cpp src/particle-model-factory.cpp src/particle-spec.cpp
        src/particle-spec-catalog.cpp src/pfactory.cpp src/polarizable-water.cpp src/protonation-site-catalog.cpp)

add_subdirectory(tests)

install(TARGETS particles DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
install(DIRECTORY resources DESTINATION ${CMAKE_INSTALL_PREFIX})
