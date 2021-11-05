cmake_minimum_required(VERSION 3.19)
project(simulation)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

include_directories(include)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)


add_library(simulation SHARED
        src/acid-base-solution.cpp src/analysis.cpp src/cell-lists.cpp src/cg-electrolyte.cpp src/cg-hp.cpp
        src/cg-lj-fluid.cpp src/cg-pol-water.cpp src/constant-rate-pt.cpp src/distance-lists.cpp
        src/interactor.cpp src/langevin-leap-frog.cpp src/langevin-velocity-verlet.cpp src/leap-frog.cpp
        src/lj-coulomb-forces.cpp src/mc.cpp src/no-bc.cpp src/pbc.cpp src/pt-langevin-velocity-verlet.cpp
        src/pt-pair-list-generator.cpp src/s-factory.cpp src/sim-data.cpp src/sim-model.cpp src/sim-model-factory.cpp
        src/sim-util.cpp src/simulation.cpp src/velocity-verlet.cpp src/protonatable-particle-model-factory.cpp
        src/continuous.cpp src/discrete.cpp)

add_subdirectory(tests)

install(TARGETS simulation DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)

