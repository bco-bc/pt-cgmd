cmake_minimum_required(VERSION 3.20)
project(bem)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

include_directories(usr/local/include)
include_directories(include)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)

add_library(bem SHARED
            src/kernels.cpp
            src/flat-tri-nodes-vertices-calculator.cpp src/flat-tri-nodes-tri-calculator.cpp
            src/curved-tri-nodes-vertices-calculator.cpp
            src/bem-data.cpp src/rhs.cpp
            src/curve.cpp)

add_subdirectory(tests)

install(DIRECTORY include/simploce DESTINATION include)
