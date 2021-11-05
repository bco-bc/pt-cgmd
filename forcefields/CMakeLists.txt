cmake_minimum_required(VERSION 3.19)
project(forcefields)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

include_directories(include)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)

add_library(forcefields SHARED
        src/no_bc.cpp)

add_subdirectory(tests)

install(TARGETS forcefields DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)

