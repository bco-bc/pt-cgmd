/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   nint-test.cpp
 * Author: juffer
 *
 * Created on 19 October 2019, 15:50
 */

#include "simploce/util/util.hpp"
#include "simploce/util/utypes.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "nint-test test 1" << std::endl;
    real_t x0 = -1.55;
    real_t dx = 0.05;
    std::size_t N = 44;
    for (std::size_t k = 0; k != N; ++k) {
        real_t x = x0 + k * dx;
        std::cout << x << ' ' << util::nint(x) << std::endl;
    }
    std::cout << std::endl;
    
    real_t a = -9.5000000e-01;
    real_t n = util::nint(a);
    std::cout << a << ' ' << n << std::endl;
}

void test2() {
    std::cout << "nint-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (nint-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% nint-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (nint-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (nint-test)" << std::endl;

//    std::cout << "%TEST_STARTED% test2 (nint-test)\n" << std::endl;
//    test2();
//    std::cout << "%TEST_FINISHED% time=0 test2 (nint-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

