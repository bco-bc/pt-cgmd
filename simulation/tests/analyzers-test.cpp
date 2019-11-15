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
 * File:   analyzers-test.cpp
 * Author: ajuffer
 *
 * Created on September 23, 2019, 3:20 PM
 */

#include "simploce/analysis/gr.hpp"
#include "simploce/analysis/dipole-moment.hpp"
#include "simploce/analysis/diffusion.hpp"
#include "simploce/analysis/analysis.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/simulation/sall.hpp"
#include "simploce/simulation/stypes.hpp"
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace simploce;

static std::vector<bead_ptr_t> all_{};
static std::vector<bead_ptr_t> free_{};
static std::vector<bead_group_ptr_t> groups_{};

void setup()
{
    CoarseGrained cg;
    auto spec = ParticleSpec::create("DW", 1.0, 10.0,  0.4);
    auto bead = cg.addBead(1, "DW", position_t{}, spec);
    all_.push_back(bead);
}

/*
 * Simple C++ Test Suite
 */

void test1() 
{
    std::cout << "analyzers-test test 1" << std::endl;
    
    using gr_t = Gr<Bead>;    
    
    length_t dr{0.1};
    std::string specName = "DW";
    box_ptr_t box = std::make_shared<box_t>(7);
    bc_ptr_t bc = factory::pbc(box);
    
    gr_t gr(dr, specName, specName, box, bc);
    gr.perform(all_, free_, groups_);
    
    cg_analyzer_ptr_t ptr = gr_t::create(dr, specName, specName, box, bc);    
}

void test2() 
{
    std::cout << "analyzers-test test 2" << std::endl;
    
    using m_m2_t = DipoleMoment<Bead>;
    
    m_m2_t mm2(10*0.02, 0.01, 100);
    mm2.perform(all_, free_, groups_);
}

void test3()
{
    std::cout << "analyzers-test test 3" << std::endl;
    using diffusion_t = Diffusion<Bead>;
    
    real_t dt = 10.0 * 0.020;
    real_t tau = 10.0;
    auto diffusion = diffusion_t::create(dt, tau);
    diffusion->perform(all_, free_, groups_);
    diffusion->results();
    
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% analyzers-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;
    
    setup();

    std::cout << "%TEST_STARTED% test1 (analyzers-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (analyzers-test)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (analyzers-test)" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (analyzers-test)" << std::endl;

    std::cout << "%TEST_STARTED% test3 (analyzers-test)" << std::endl;
    test3();
    std::cout << "%TEST_FINISHED% time=0 test3 (analyzers-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

