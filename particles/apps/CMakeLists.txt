add_executable(total-masses total-masses.cpp)

target_link_directories(total-masses PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(total-masses particles util cpputil)

install(TARGETS total-masses DESTINATION bin)