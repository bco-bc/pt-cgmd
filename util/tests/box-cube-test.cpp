/*
 * File:   box-cube-test.cpp
 * Author: André Juffer
 *
 * Created on August 29, 2019, 5:48 PM
 */

#include "simploce/util/cube.hpp"
#include <iostream>

using namespace simploce;

int main(int argc, char** argv) {
    Box<double> box{1.0, 2.0, 3.0};
    Cube<double> cube{3.5};
    std::cout << "Box. " << std::endl << box << std::endl;
    std::cout << "Cube: " << std::endl << cube << std::endl;
    std::cout << "Cube: " << std::endl;
    for (auto k = 0; k !=3; ++k) {
        std::cout << cube[k] << std::endl;
    }
    std::cout << "Box size: " << box.size() << std::endl;
    std::cout << "Cube size: " << cube.size() << std::endl;

    std::cout << "Read cuve (box): " << std::endl;
    std::cin >> cube;
    std::cout << "Cube: " << cube << std::endl;
    return EXIT_SUCCESS;
}

