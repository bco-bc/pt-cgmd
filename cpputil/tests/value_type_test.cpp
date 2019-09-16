/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   value_type_test.cpp
 * Author: ajuffer
 *
 * Created on August 9, 2019, 1:37 PM
 */

#include "simploce/util/value_t.hpp"
#include <cstdlib>
#include <iostream>

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "value_type_test test 1" << std::endl;
    
    using result_t = simploce::value_t<double, 0>;
    
    using pressure_t = simploce::value_t<double, 1>;
    using volume_t = simploce::value_t<double, 2>;
    using temperature_t = simploce::value_t<double, 3>;
    
    pressure_t p1{2};
    pressure_t p3 = p1;
    pressure_t p2(4);
    pressure_t p4{p2};
    std::cout << "p1: " << p1 << std::endl;
    std::cout << "p2: " << p2 << std::endl;
    std::cout << "p3 = p1: " << p3 << std::endl;
    std::cout << "p4{p2}: " << p4 << std::endl;
    std::cout << "(p2 == p4)? " << (p2==p4) << std::endl;
    std::cout << "(p2 == p3)? " << (p2==p3) << std::endl;
    pressure_t p = p1 + p2;
    std::cout << "p : " << p << std::endl;
        
    volume_t V(300);
    std::cout << "V : " << V << std::endl;
    pressure_t r1 = p * 2.0;
    std::cout << "p * 2.0: " << r1 << std::endl;    
    std::cout << "3.0 * V: " << 3.0 * V << std::endl;
    
    result_t r2 = p * V;
    std::cout << "p * V: " << r2 << std::endl;
    std::cout << "p()*V(): " << p() * V() << std::endl;
    
    std::cout << "2.0 * V * 3.0  : " << 2.0 * V * 3.0 << std::endl;
    std::cout << "2.0 * V() * 3.0: " << 2.0 * V() * 3.0 << std::endl;
    
    temperature_t T(300);
    
    std::cout << "p * V / T: " << p * V / T << std::endl;
    std::cout << "p() * V() / T(): " << p() * V() / T() << std::endl;
    
    // Must fail.
    //result_t r4{p};
    //std::cout << "p + V: " << p + V << std::endl;
}

void test2() {
    std::cout << "value_type_test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (value_type_test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% value_type_test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (value_type_test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (value_type_test)" << std::endl;

    /*
    std::cout << "%TEST_STARTED% test2 (value_type_test)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (value_type_test)" << std::endl;
    */
    
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

