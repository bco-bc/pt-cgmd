cmake_minimum_required(VERSION 3.19)
project(apps)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX /localdisk)
if(PREFIX)
    message("PREFIX=${PREFIX}")
    set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

# Threading
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Enable dynamic linking with boost libraries.
set(Boost_USE_STATIC_LIBS OFF)

# Required boost library
find_package(Boost REQUIRED COMPONENTS program_options)

include_directories(${Boost_INCLUDE_DIR} SYSTEM)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)

#[[
add_executable(s-dipole-moment src/s-dipole-moment.cpp)
add_executable(s-msd src/s-msd.cpp)]]
add_executable(s-particle-system src/s-particle-system.cpp)
add_executable(s-validate src/s-validate.cpp)
add_executable(s-simulate src/s-simulate.cpp)
add_executable(s-gr src/s-gr.cpp)
add_executable(s-distance src/s-distance.cpp)

#[[target_link_directories(s-dipole-moment PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-dipole-moment simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-msd PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-msd simulation particles util cpputil pthread ${Boost_LIBRARIES})]]
target_link_directories(s-particle-system PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-particle-system simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-validate PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-validate simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-simulate PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-simulate simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-gr PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-gr simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-distance PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-distance simulation particles util cpputil pthread ${Boost_LIBRARIES})

#[[install(TARGETS s-cgmd s-dipole-moment s-gr s-msd s-particle-system DESTINATION bin)]]
install(TARGETS s-particle-system s-validate s-simulate s-gr s-distance DESTINATION bin)