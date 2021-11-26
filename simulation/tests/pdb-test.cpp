/*
 * File:   pdb-test.cpp
 * Author: ajuffer
 *
 * Created on October 18, 2019, 3:13 PM
 */

#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/types/cvector_t.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "pdb-test test 1" << std::endl;
    
    box_ptr_t box = factory::box(2.0);
    std::cout << "Box size: " << box->size() << std::endl;
    bc_ptr_t bc = factory::boundaryCondition(box);
    
    position_t r1{0,0,0}; 
    real_t dz = 0.1;
    std::size_t N = 20;
    
    std::cout.setf(std::ios::scientific);
    for (std::size_t k = 1; k != N; ++k) {
        real_t z = k * dz;
        position_t r2{0, 0, z};
        length_t R = norm<real_t>(r1-r2);
        auto r12 = bc->apply(r1, r2);
        length_t Rbc = norm<real_t>(r12);
        std::cout << "k, R, Rbc: " << k << ' ' << R << ' ' << Rbc << std::endl;
    }
}

void test2() {
    std::cout << "pdb-test test 2" << std::endl;
    
    box_ptr_t box = factory::box(2.0);
    std::cout << "Box size: " << box->size() << std::endl;
    bc_ptr_t bc = factory::boundaryCondition(box);
    
    real_t z0 = -0.5;
    real_t dz = 0.1;
    std::size_t N = 31;
    
    std::cout.setf(std::ios::scientific);
    for (std::size_t k = 0; k != N; ++k) {
        real_t z = z0 + k * dz;
        position_t r{z, z, z};
        auto rin = bc->placeInside(r);
        std::clog << "Actual: " << r << ", Placed: " << rin << std::endl;
    }
}


int main(int argc, char** argv) {
    test1();
    test2();
    return (EXIT_SUCCESS);
}

