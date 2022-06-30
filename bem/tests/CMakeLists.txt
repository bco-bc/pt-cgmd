set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(test-matrix-vector test-matrix-vector.cpp)
add_executable(test-kernels test-kernels.cpp)
add_executable(test-bem-calculator test-bem-calculator.cpp)

target_link_directories(test-matrix-vector PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(test-matrix-vector bem surface util cpputil pthread)
target_link_directories(test-kernels PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(test-kernels bem surface util cpputil pthread)
target_link_directories(test-bem-calculator PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(test-bem-calculator bem surface util cpputil pthread)
