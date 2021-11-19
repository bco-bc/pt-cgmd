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

#[[add_library(simulation SHARED
        src/acid-base-solution.cpp src/analysis.cpp src/cell-lists.cpp src/cg-electrolyte.cpp src/cg-hp.cpp
        src/cg-lj-fluid.cpp src/cg-pol-water.cpp src/constant-rate-pt.cpp src/distance-lists.cpp
        src/interactor.cpp src/langevin-leap-frog.cpp src/langevin-velocity-verlet.cpp src/leap-frog.cpp
        src/lj-coulomb-forces.cpp src/mc.cpp src/no-bc.cpp src/pbc.cpp src/pt-langevin-velocity-verlet.cpp
        src/pt-pair-list-generator.cpp src/s-factory.cpp src/sim-data.cpp src/sim-model.cpp src/sim-model-factory.cpp
        src/s-properties.cpp src/simulation.cpp src/velocity-verlet.cpp src/protonatable-particle-system-factory.cpp
        src/continuous.cpp src/discrete.cpp src/force-field.cpp src/distance-pair-list-generator.cpp)
]]

add_library(simulation SHARED
        src/interactor.cpp src/forces.cpp src/lj.cpp src/lj-rf.cpp src/hp.cpp
        src/halve-attractive-qp.cpp src/pair-lists.cpp src/simulation.cpp src/sim-data.cpp
        src/no-bc.cpp src/pbc.cpp src/s-factory.cpp
        src/s-properties.cpp src/protonatable-particle-system-factory.cpp
        src/continuous.cpp src/discrete.cpp src/force-field.cpp
        src/distance-pair-list-generator.cpp src/leap-frog.cpp src/velocity-verlet.cpp
        src/langevin-leap-frog.cpp src/langevin-velocity-verlet.cpp
        src/mc.cpp)


add_subdirectory(tests)
add_subdirectory(apps)

install(TARGETS simulation DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
install(DIRECTORY resources DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS lj-parameters DESTINATION bin)
