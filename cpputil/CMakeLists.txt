cmake_minimum_required(VERSION 3.19)
project(cpputil)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /home/andre)

include_directories(include)
add_library(cpputil SHARED src/file.cpp src/param.cpp src/poisson-process.cpp src/telegraph-process.cpp src/util.cpp)
add_subdirectory(tests)

install(TARGETS cpputil DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
