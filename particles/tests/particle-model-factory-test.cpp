/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particle-model-factory-test.cpp
 * Author: ajuffer
 *
 * Created on October 3, 2019, 1:31 PM
 */

#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/pfactory.hpp"
#include "simploce/particle/ptypes.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "particle-model-factory-test test 1" << std::endl;
    
    std::string fileName = 
            "/home/ajuffer/simploce/pt-cgmd/particles/resources/particles-specs.dat";
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fileName);
    std::clog << *catalog << std::endl;
    
    particle_model_fact_ptr_t factory = factory::particleModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{3.0});
    cg_ptr_t cg = factory->formicAcidSolution(box);
    std::clog << *cg << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% particle-model-factory-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (particle-model-factory-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (particle-model-factory-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

