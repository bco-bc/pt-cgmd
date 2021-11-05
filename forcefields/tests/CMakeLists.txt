add_executable(ext-potential-test ext-potential-test.cpp)
add_executable(static-uniform-test static-uniform-test.cpp)
add_executable(morse-test morse-test.cpp)
add_executable(pair-potentials-test pair-potentials-test.cpp)

target_link_directories(ext-potential-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(ext-potential-test forcefields particles util cpputil)
target_link_directories(static-uniform-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(static-uniform-test forcefields particles util cpputil)
target_link_directories(morse-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(morse-test forcefields particles util cpputil)
target_link_directories(pair-potentials-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(pair-potentials-test forcefields particles util cpputil)

