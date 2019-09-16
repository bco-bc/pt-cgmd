/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   value-cvector-test.cpp
 * Author: ajuffer
 *
 * Created on August 16, 2019, 1:12 PM
 */

#include <cstdlib>
#include <iostream>

#include "simploce/util/cvector_t.hpp"
#include "simploce/util/value_t.hpp"

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "value-cvector-test test 1" << std::endl;
    
    using momentum_t = simploce::cvector_t<double, 1>;
    using velocity_t = simploce::cvector_t<double, 2>;
    using position_t = simploce::cvector_t<double, 3>;
    using force_t = simploce::cvector_t<double, 4>;
    using stime_t = simploce::value_t<double, 1>;
    using mass_t = simploce::value_t<double, 2>;
    
    velocity_t vi, vf;
    momentum_t p{1,2,3};
    mass_t mass{2.0};
    stime_t dt = 0.001;
    force_t f{9,7,5};
    position_t r;
    
    for (std::size_t k = 0; k != 3; ++k) {
        vi[k] = p[k] / mass();                     // Velocity (nm/ps) at time t(n-1/2)
        p[k] += dt() * f[k];                      // Momentum at time t(n+1/2)
        
        vf[k] = p[k] / mass();                     // Velocity at time t(n+1/2)
        r[k] += dt() * vf[k];                     // Position at time t(n+1).
    }
    
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% value-cvector-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (value-cvector-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (value-cvector-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

