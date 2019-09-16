/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   seed-value-test.cpp
 * Author: ajuffer
 *
 * Created on August 21, 2019, 1:37 PM
 */

#include "simploce/util/util.hpp"
#include <cstdlib>
#include <iostream>

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "seed-value-test test 1" << std::endl;
    
    for ( std::size_t counter = 0; counter != 10; ++counter) {
        auto value = simploce::util::seedValue<std::size_t>();
        std::cout << "Seed value: " << value << std::endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% seed-value-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (seed-value-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (seed-value-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

