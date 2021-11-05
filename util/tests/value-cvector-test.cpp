/*
 * File:   value-cvector-test.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 1:12 PM
 */

#include <cstdlib>
#include <iostream>

#include "simploce/types/u-types.hpp"

int main(int argc, char** argv) {

    simploce::velocity_t vi, vf;
    simploce::momentum_t p{1,2,3};
    simploce::mass_t mass{2.0};
    simploce::stime_t dt = 0.001;
    simploce::force_t f{9000,7000,5000};
    simploce::position_t ri, rf;
    simploce::energy_t ekin{100.0}, epot(-212);

    std::cout << "Before:" << std::endl;
    std::cout << "Velocity: " << vi << std::endl;
    for (std::size_t k = 0; k != 3; ++k) {
        vi[k] = p[k] / mass();                     // Velocity (nm/ps) at time t(n-1/2)
        p[k] += dt() * f[k];                       // Momentum at time t(n+1/2)

        vf[k] = p[k] / mass();                     // Velocity at time t(n+1/2)
        rf[k] += ri[k] + dt() * vf[k];             // Position at time t(n+1).
    }
    std::cout << "Initial:" << std::endl;
    std::cout << "Velocity: " << vi << std::endl;
    std::cout << "Position: " << ri << std::endl;
    std::cout << "Final:" << std::endl;
    std::cout << "Velocity: " << vf << std::endl;
    std::cout << "Position: " << rf << std::endl;

    std::cout << "Energy: " << 0.5 * (ekin + epot) << std::endl;
    return (EXIT_SUCCESS);
}

