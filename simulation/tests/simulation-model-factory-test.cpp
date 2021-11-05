/*
 * File:   model-factory-test.cpp
 * Author: ajuffer
 *
 * Created on August 30, 2019, 1:53 PM
 */

#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-types.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/cube.hpp"
#include "simploce/simulation/s-conf.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() 
{
    std::cout << "simulation-model-factory-test test 1" << std::endl;
    std::string fileName = "resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
    std::cout << *catalog << std::endl;
    
    sim_model_fact_ptr_t pmf = 
            factory::simulationModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{7.27});
    cg_sim_model_ptr_t polWater = 
        pmf->polarizableWater(box, 997.0479, 298.15, 2);

    std::cout << "Number of polarizable waters: " << polWater->size() << std::endl;
    std::cout << *polWater << std::endl;
    
}

void test2() 
{
    std::cout << "simulation-model-factory-test test 2" << std::endl;
    
    std::string fileName = "resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
    //std::cout << *catalog << std::endl;
    
    sim_model_fact_ptr_t pmf = 
            factory::simulationModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{7.27});
    cg_sim_model_ptr_t hcooh = pmf->formicAcidSolution(box);
    
    std::cout << "Number of beads in HCCOH model: " << hcooh->size() << std::endl;
    std::cout << *hcooh << std::endl;
}

void test3()
{
    std::cout << "simulation-model-factory-test test 2" << std::endl;
    
    std::string fileName = "resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
    //std::cout << *catalog << std::endl;
    
    sim_model_fact_ptr_t pmf = 
            factory::simulationModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{7.27});
    cg_sim_model_ptr_t nacl = pmf->electrolyte(box);
    
    std::cout << "Number of beads in NaCl model: " << nacl->size() << std::endl;
    std::cout << *nacl << std::endl;
}

void test4() 
{
    std::cout << "simulation-model-factory-test test 4" << std::endl;
    
    std::string fileName = "resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
    //std::cout << *catalog << std::endl;
    
    sim_model_fact_ptr_t pmf = 
            factory::simulationModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{7.27});
    cg_sim_model_ptr_t hcooh = pmf->formicAcidSolution(box);
    std::string displacerId = conf::LANGEVIN_VELOCITY_VERLET;
    
    std::cout << "Number of beads in HCCOH model: " << hcooh->size() << std::endl;
    std::cout << *hcooh << std::endl;
}

void test5() 
{
    std::cout << "simulation-model-factory-test test 5" << std::endl;
    
    std::string fileName = "resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
    //std::cout << *catalog << std::endl;
    
    sim_model_fact_ptr_t pmf = 
            factory::simulationModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{7.27});
    cg_sim_model_ptr_t sm = pmf->ljFluid(box);
    std::string displacerId = conf::LANGEVIN_VELOCITY_VERLET;
    
    std::cout << "Number of beads in LJ fluid model: " << sm->size() << std::endl;
    std::cout << *sm << std::endl;
}



int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% model-factory-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (model-factory-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (model-factory-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test2 (model-factory-test)\n" << std::endl;
    //test2();
    //std::cout << "%TEST_FINISHED% time=0 test2 (model-factory-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test3 (model-factory-test)\n" << std::endl;
    //test3();
    //std::cout << "%TEST_FINISHED% time=0 test3 (model-factory-test)" << std::endl;
    
    //std::cout << "%TEST_STARTED% test4 (model-factory-test)\n" << std::endl;
    //test4();
    //std::cout << "%TEST_FINISHED% time=0 test4 (model-factory-test)" << std::endl;
    
    std::cout << "%TEST_STARTED% test5 (model-factory-test)\n" << std::endl;
    test5();
    std::cout << "%TEST_FINISHED% time=0 test5 (model-factory-test)" << std::endl;
    
    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

