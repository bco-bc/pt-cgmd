cmake_minimum_required(VERSION 3.20)
project(dw)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /home/juffer)

include_directories(/usr/local/include)
include_directories(include)

install(DIRECTORY include/simploce DESTINATION include)
