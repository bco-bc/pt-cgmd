/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   pdb-test.cpp
 * Author: ajuffer
 *
 * Created on October 18, 2019, 3:13 PM
 */

#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/util/cvector_t.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "pdb-test test 1" << std::endl;
    
    box_ptr_t box = factory::cube(2.0);
    std::cout << "Box size: " << box->size() << std::endl;
    bc_ptr_t bc = factory::pbc(box);
    
    position_t r1{0,0,0}; 
    real_t dz = 0.1;
    std::size_t N = 20;
    
    for (std::size_t k = 1; k != N; ++k) {
        real_t z = k * dz;
        position_t r2{0, 0, z};
        length_t R = norm<real_t>(r1-r2);
        auto r12 = bc->apply(r1, r2);
        length_t Rbc = norm<real_t>(r12);
        std::cout << "k, R, Rbc: " << k << ' ' << R << ' ' << Rbc << std::endl;
    }
}


int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% pdb-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (pdb-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (pdb-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

