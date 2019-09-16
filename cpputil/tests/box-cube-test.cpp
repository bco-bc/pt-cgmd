/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   box-cube-test.cpp
 * Author: ajuffer
 *
 * Created on August 29, 2019, 5:48 PM
 */

#include "simploce/util/cube.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "box-cube-test test 1" << std::endl;
    
    Box<double> box{1.0, 2.0, 3.0};
    Cube<double> cube{3.5};
    std::cout << box << std::endl;
    std::cout << cube << std::endl;
    for (auto k = 0; k !=3; ++k) {
        std::cout << cube[k] << std::endl;
    }
    std::cout << "Box size: " << box.size() << std::endl;
    std::cout << "Cube size: " << cube.size() << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% box-cube-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (box-cube-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (box-cube-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

