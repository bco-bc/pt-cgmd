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
 * File:   gr-test.cpp
 * Author: ajuffer
 *
 * Created on September 23, 2019, 3:20 PM
 */

#include "simploce/analysis/gr.hpp"
#include "simploce/analysis/analysis.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/simulation/sall.hpp"
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "gr-test test 1" << std::endl;
    
    using gr_t = Gr<Bead>;    
    
    length_t dr{0.1};
    std::string specName = "DW";
    box_ptr_t box = std::make_shared<box_t>(7);
    bc_ptr_t bc = factory::pbc(box);
    
    gr_t gr(dr, specName, specName, box, bc);
    
    cg_analyzer_ptr_t ptr = gr_t::create(dr, specName, specName, box, bc);    
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% gr-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (gr-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (gr-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

