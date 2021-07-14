/*
 * File:   simulation-test.cpp
 * Author: ajuffer
 *
 * Created on August 20, 2019, 12:50 PM
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/simulation/cg-forcefield.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>

using namespace simploce;

struct SimpleForceField : public CoarseGrainedForceField{
    
    std::pair<energy_t, energy_t>
    interact(const std::vector<bead_ptr_t>& all,
             const std::vector<bead_ptr_t>& free,
             const std::vector<bead_group_ptr_t>& groups,
             const PairLists<Bead>& pairLists) override    
    {
        std::clog << "Computing forces for simulation for all, free, and groups...Done" << std::endl;
        return std::make_pair(0.0, 0.0);       
    }
    
    energy_t bonded(const std::vector<bead_ptr_t>& all,
                    const std::vector<bead_ptr_t>& free,
                    const std::vector<bead_group_ptr_t>& groups,
                    const PairLists<Bead>& pairLists)  override {
        return 0.0;
    }
    
    std::pair<energy_t, energy_t>
    interact(const bead_ptr_t& bead,
             const std::vector<bead_ptr_t>& all,
             const std::vector<bead_ptr_t>& free,
             const std::vector<bead_group_ptr_t>& groups) override
    {
        return std::make_pair(0.0, 0.0);     
    }
    
    std::string id() const { return "cg-simple-ff"; }
    
    std::pair<lj_params_t, el_params_t> parameters() const override { 
        return std::pair<lj_params_t, el_params_t>{}; 
    }
};


/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "simulation-test test 1" << std::endl;
    
    std::ofstream trajStream, dataStream;
    file::open_output(trajStream, "tests/trajectory.dat");
    file::open_output(dataStream, "tests/simulation.dat");
    
    using simulation_t = Simulation<Bead>;
    //using simulation_t = MC<Bead>;

    spec_ptr_t spec1_dp = 
        ParticleSpec::create("pspec1", 1.0, 1.0, 1.0);
    spec_ptr_t spec1_p = ParticleSpec::create("dpspec1", 1.0, 1.0, 1.0);
    
    cg_ptr_t cg = factory::coarseGrained();
    cg->addBead(1, "test1", position_t{}, spec1_dp);
    cg->addBead(2, "test2", position_t{}, spec1_p);
    std::cout << "Number of beads: " << cg->numberOfParticles() << std::endl;
    
    box_ptr_t box = factory::cube(length_t{1.0});
    bc_ptr_t bc = factory::pbc(box);
    
    cg_ff_ptr_t forcefield = std::make_shared<SimpleForceField>();
    cg_ppair_list_gen_ptr_t generator = factory::coarseGrainedPairListGenerator(box, bc);
    cg_interactor_ptr_t interactor = std::make_shared<Interactor<Bead>>(forcefield, generator);
    pt_displacer_ptr_t ptDisplacer = factory::protonTransferDisplacer();
    pt_pair_list_gen_ptr_t ptGenerator = factory::protonTransferPairListGenerator(bc);
    cg_displacer_ptr_t displacer = 
            factory::protonTransferlangevinVelocityVerlet(interactor, 
                                                          ptGenerator, 
                                                          ptDisplacer);
    
    cg_sim_model_ptr_t sm = 
        std::make_shared<cg_sim_model_t>(cg, displacer, interactor, box, bc);
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

