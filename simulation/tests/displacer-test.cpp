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
 * File:   leap-frog-test.cpp
 * Author: ajuffer
 *
 * Created on August 16, 2019, 1:00 PM
 */

#include <cstdlib>
#include <iostream>


#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/stypes.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"

using namespace simploce;


/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "displacer-test test 1" << std::endl;
    
    sim_param_t param{};
    param.add<std::size_t>("nsteps", 1000);
    param.add<real_t>("timestep", 0.001);
    
    cg_ptr_t cg = factory::coarseGrained();

    box_ptr_t box = factory::cube(length_t{5.0});
    bc_ptr_t bc = factory::pbc(box);
    cg_interactor_ptr_t interactor = 
        factory::interactorCoarseGrainedPolarizableWater(box, bc);
    
    cg_displacer_ptr_t leapFrog = factory::leapFrog(interactor);
    SimulationData data = leapFrog->displace(param, cg);
    
    cg_displacer_ptr_t velocityVerlet = factory::velocityVerlet(interactor);
    data = velocityVerlet->displace(param, cg);
    std::cout << data << std::endl;
          
    cg_displacer_ptr_t lvelocityVerlet = factory::langevinVelocityVerlet(interactor);
    data = lvelocityVerlet->displace(param, cg);
    std::cout << data << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% displacer-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (displacer-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (displacer-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;
    
    return (EXIT_SUCCESS);
}

