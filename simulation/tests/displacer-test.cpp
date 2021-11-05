/*
 * File:   leap-frog-test.cpp
 * Author: ajuffer
 *
 * Created on August 16, 2019, 1:00 PM
 */

#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-types.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/param.hpp"
#include "simploce/simulation/sim-model-factory.hpp"
#include <fstream>
#include <cstdlib>
#include <iostream>

using namespace simploce;
using namespace simploce::param;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "displacer-test test 1" << std::endl;
    
    std::string fileName = "resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);

    sim_param_t param{};
    param.add<std::size_t>("nsteps", 1000);
    param.add<real_t>("timestep", 0.001);
    param.add<real_t>("npairlists", 10);
    param.add<real_t>("temperature", 300);
    param.add<real_t>("gamma", 0.5);
    std::cout << param << std::endl;
    
    box_ptr_t box = factory::cube(length_t{5.0});
    bc_ptr_t bc = factory::pbc(box);

    particle_model_fact_ptr_t factory = factory::particleModelFactory(catalog);
    CoarseGrained cg = factory->polarizableWater(box);
    cg_interactor_ptr_t interactor = 
            factory::polarizableWaterInteractor(catalog, box, bc);
    
    /*
    cg_displacer_ptr_t leapFrog = factory::leapFrog(interactor);    
    SimulationData data = leapFrog->displace(param, cg);
    
    cg_displacer_ptr_t velocityVerlet = factory::velocityVerlet(interactor);
    
    data = velocityVerlet->displace(param, cg);
    std::cout << data << std::endl;
          
    cg_displacer_ptr_t lvelocityVerlet = factory::langevinVelocityVerlet(interactor);
    data = lvelocityVerlet->displace(param, cg);
    std::cout << data << std::endl;
    */
    /*
    cg_ptr_t ptcg = factory->formicAcidSolution(box);
    pt_pair_list_gen_ptr_t generator = factory::protonTransferPairListGenerator(bc);
    pt_displacer_ptr_t ptDisplacer = factory::protonTransferDisplacer();
    cg_displacer_ptr_t pt = factory::protonTransferlangevinVelocityVerlet(interactor,
                                                                          generator,
                                                                          ptDisplacer);
    SimulationData data = pt->displace(param, ptcg);
    std::cout << data << std::endl;
     */
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

