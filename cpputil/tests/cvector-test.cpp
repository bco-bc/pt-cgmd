/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cvector-test.cpp
 * Author: ajuffer
 *
 * Created on August 8, 2019, 12:03 PM
 */

#include "simploce/util/cvector_t.hpp"
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "cvector-test test 1" << std::endl;
    
    using iv_t = simploce::cvector_t<int, 2>;
   
    iv_t iv{1, 2, 3};
    std::cout << "iv: " << iv << std::endl;
    std::cout << "|iv|: " << simploce::norm<double>(iv) << std::endl;
    
    std::array<double, 3> av{10.0, 11.0, -12.0};
    
    using position_t = simploce::cvector_t<double, 1>;
    position_t r1{1.0, 0, 0};
    position_t r2{0.0, 0.0, 2.0};
    position_t r3{av};
    
    std::cout << "r1: " << r1 << std::endl;
    std::cout << "r2: " << r2 << std::endl;
    std::cout << "r2: " << r3 << std::endl;
    
    position_t rd = r2 - r1 + r3;
    std::cout << "rd: " << rd << std::endl;
    std::cout << "|rd| = " << simploce::norm<double>(rd) << std::endl;
    
    rd *= 4.0;
    std::cout << "rd: " << rd << std::endl;
    
    auto a = simploce::angle<double>(r1, r2);
    std::cout << "Angle (r1, r2): " << a << std::endl;
    std::cout << "0.5 * PI      : " << 0.5 * std::acos(-1.L) << std::endl;
    
    auto cp = simploce::cross<double,2>(r1, r2);
    std::cout << "cp: " << cp << std::endl;
    std::cout << "Angle (r1, cp): " << simploce::angle<double>(r1, cp) << std::endl;
}

void test2() {
    std::cout << "cvector-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (cvector-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% cvector-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (cvector-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (cvector-test)" << std::endl;

    /*
    std::cout << "%TEST_STARTED% test2 (cvector-test)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (cvector-test)" << std::endl;
    */
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

