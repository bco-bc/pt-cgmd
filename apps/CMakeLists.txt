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

include_directories(include)
include_directories(${Boost_INCLUDE_DIR} SYSTEM)
include_directories(${CMAKE_INSTALL_PREFIX}/include SYSTEM)

add_library(apps SHARED src/program-options.cpp)

add_executable(s-validate src/s-validate.cpp)
add_executable(s-simulate src/s-simulate.cpp)
add_executable(s-gr src/s-gr.cpp)
add_executable(s-particle-system src/s-particle-system.cpp)
add_executable(s-pdb src/s-pdb.cpp)
add_executable(s-add-boundary src/s-add-boundary.cpp)
add_executable(s-triangulate src/s-triangulate.cpp)

target_link_directories(s-validate PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-validate apps simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-simulate PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-simulate apps simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-gr PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-gr apps simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-particle-system PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-particle-system apps simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-pdb PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-pdb apps simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-add-boundary PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-add-boundary apps simulation particles util cpputil pthread ${Boost_LIBRARIES})
target_link_directories(s-triangulate PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(s-triangulate apps surface simulation particles util cpputil pthread ${Boost_LIBRARIES})

install(TARGETS apps DESTINATION lib)
install(DIRECTORY include/simploce DESTINATION include)
install(TARGETS s-particle-system s-validate s-simulate s-gr s-pdb s-add-boundary s-triangulate DESTINATION bin)