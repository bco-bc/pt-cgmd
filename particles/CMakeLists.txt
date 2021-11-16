cmake_minimum_required(VERSION 3.19)
project(particles)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

include_directories(include)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)

add_library(particles SHARED
        src/atom.cpp src/atomistic.cpp src/bead.cpp src/coarse-grained.cpp
        src/particle.cpp src/particle-system-factory.cpp src/particle-spec.cpp
        src/particle-spec-catalog.cpp src/p-factory.cpp
        src/protonation-site-catalog.cpp)

add_subdirectory(apps)
add_subdirectory(tests)

install(TARGETS particles DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
install(DIRECTORY resources DESTINATION ${CMAKE_INSTALL_PREFIX})
