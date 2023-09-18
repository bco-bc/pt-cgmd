add_executable(lj-parameters lj-parameters.cpp)
add_executable(one-water-dpd-particle-simulation one-water-dpd-particle-simulation.cpp)

target_link_directories(lj-parameters PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(lj-parameters simulation particles util cpputil pthread)
target_link_directories(one-water-dpd-particle-simulation PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(one-water-dpd-particle-simulation simulation particles util cpputil pthread)
