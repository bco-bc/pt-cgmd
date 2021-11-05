cmake_minimum_required(VERSION 3.19)
project(util)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

include_directories(include)
add_library(util SHARED
        src/param.cpp src/poisson-process.cpp src/telegraph-process.cpp
        src/util.cpp)
add_subdirectory(tests)

install(TARGETS util DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
