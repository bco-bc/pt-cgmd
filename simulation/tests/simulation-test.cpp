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
 * File:   simulation-test.cpp
 * Author: ajuffer
 *
 * Created on August 20, 2019, 12:50 PM
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/simulation/cg-forcefield.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>

using namespace simploce;

struct SimpleForceField : public CoarseGrainedForceField{
    
    energy_t interact(const std::vector<bead_ptr_t>& all,
                      const std::vector<bead_ptr_t>& free,
                      const std::vector<bead_group_ptr_t>& groups,
                      const std::vector<bead_pair_list_t>& pairLists) override    
    {
        std::clog << "Computing forces for simulation for all, free, and groups...Done" << std::endl;
        return 0.0;        
    }
    
    std::string id() const { return "cg-simple-ff"; }
};


/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "simulation-test test 1" << std::endl;
    
    std::ofstream trajStream, dataStream;
    file::open_output(trajStream, "/localdisk/tests/trajectory.dat");
    file::open_output(dataStream, "/localdisk/tests/simulation.dat");
    
    using simulation_t = Simulation<Bead>;

    spec_ptr_t spec1_dp = 
        ParticleSpec::create("pspec1", 1.0, 1.0, 1.0);
    spec_ptr_t spec1_p = ParticleSpec::create("dpspec1", 1.0, 1.0, 1.0);
    
    cg_ptr_t cg = factory::coarseGrained();
    cg->addBead(1, "test1", position_t{}, spec1_dp);
    cg->addBead(2, "test2", position_t{}, spec1_p);
    std::cout << "Number of beads: " << cg->numberOfParticles() << std::endl;
    
    box_ptr_t box = factory::cube(length_t{5.0});
    bc_ptr_t bc = factory::pbc(box);
    
    cg_ff_ptr_t forcefield = std::make_shared<SimpleForceField>();
    cg_ppair_list_gen_ptr_t generator = factory::coarseGrainedPairListGenerator(box, bc);
    cg_interactor_ptr_t interactor = std::make_shared<Interactor<Bead>>(forcefield, generator);
    cg_displacer_ptr_t lvv = factory::langevinVelocityVerlet(interactor);
    
    cg_sim_model_ptr_t sm = std::make_shared<cg_sim_model_t>(cg, lvv, interactor, box, bc);
    simulation_t simulation{sm};
    sim_param_t param;
    param.add<std::size_t>("nsteps", 100);
    param.add<real_t>("timestep", 0.001);
    param.add<real_t>("temperature", 298);
    param.add<real_t>("gamma", 0.5);
    param.add<std::size_t>("nwrite", 10);
    param.add<std::size_t>("npairlists", 10);
    simulation.perform(param, trajStream, dataStream);
    
    trajStream.close();
    dataStream.close();
    
}

void test2() {
    std::cout << "simulation-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (simulation-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% simulation-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (simulation-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (simulation-test)" << std::endl;

//    std::cout << "%TEST_STARTED% test2 (simulation-test)\n" << std::endl;
//    test2();
//    std::cout << "%TEST_FINISHED% time=0 test2 (simulation-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

