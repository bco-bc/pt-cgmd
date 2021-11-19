add_executable(lj-parameters lj-parameters.cpp)

target_link_directories(lj-parameters PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(lj-parameters simulation particles util cpputil pthread)
