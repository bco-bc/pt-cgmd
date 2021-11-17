add_executable(coarse-grained-test coarse-grained-test.cpp)
add_executable(particle-system-factory-test particle-system-factory-test.cpp)
add_executable(particle-spec-catalog-test particle-spec-catalog-test.cpp)
add_executable(particle-test particle-test.cpp)
#[[add_executable(prot-site-catalog-test prot-site-catalog-test.cpp)]]
#[[add_executable(protonation-site-test protonation-site-test.cpp)]]
add_executable(read-cg-test read-cg-test.cpp)
add_executable(protonatable-bead-test protonatable-bead-test.cpp)
add_executable(protonatable-coarse-grained-test protonatable-coarse-grained-test.cpp)

target_link_directories(coarse-grained-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(coarse-grained-test particles util cpputil)
target_link_directories(particle-system-factory-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(particle-system-factory-test particles util cpputil)
target_link_directories(particle-spec-catalog-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(particle-spec-catalog-test particles util cpputil)
target_link_directories(particle-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(particle-test particles util cpputil)
#[[target_link_directories(prot-site-catalog-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(prot-site-catalog-test particles util cpputil)]]
#[[target_link_directories(protonation-site-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(protonation-site-test particles util cpputil)]]
target_link_directories(read-cg-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(read-cg-test particles util cpputil)
target_link_directories(protonatable-bead-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(protonatable-bead-test particles util cpputil)
target_link_directories(protonatable-coarse-grained-test PUBLIC ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(protonatable-coarse-grained-test particles util cpputil)
