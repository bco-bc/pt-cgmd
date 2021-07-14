cmake_minimum_required(VERSION 3.19)
project(apps)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /home/andre)

# Threading
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Enable dynamic linking with boost libraries.
set(Boost_USE_STATIC_LIBS OFF)

# Required boost library
find_package(Boost REQUIRED COMPONENTS program_options)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
include_directories(${Boost_INCLUDE_DIR})

add_executable(s-cgmd src/s-cgmd.cpp)
add_executable(s-dipole-moment src/s-dipole-moment.cpp)
add_executable(s-gr src/s-gr.cpp)
add_executable(s-msd src/s-msd.cpp)

target_link_directories(s-cgmd PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-cgmd simulation particles cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-dipole-moment PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-dipole-moment simulation particles cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-gr PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-gr simulation particles cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-msd PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-msd simulation particles cpputil pthread ${Boost_LIBRARIES})

install(TARGETS s-cgmd s-dipole-moment s-gr s-msd DESTINATION bin)