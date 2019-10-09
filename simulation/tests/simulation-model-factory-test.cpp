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
 * File:   model-factory-test.cpp
 * Author: ajuffer
 *
 * Created on August 30, 2019, 1:53 PM
 */

#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/stypes.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/cube.hpp"
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
    std::string fileName = "/home/ajuffer/simploce/pt-cgmd/particles/resources/particles-specs.dat";
    std::ifstream stream;
    file::open_input(stream, fileName);
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
    std::cout << *catalog << std::endl;
    
    sim_model_fact_ptr_t pmf = 
            factory::simulationModelFactory(catalog);
    
    box_ptr_t box = factory::cube(length_t{7.27});
    cg_sim_model_ptr_t polWater = pmf->polarizableWater(box);

    std::cout << "Number of polarizable waters: " << polWater->size() << std::endl;
    //std::cout << *polWater << std::endl;
    
}

void test2() 
{
    std::cout << "simulation-model-factory-test test 2" << std::endl;
    
    std::string fileName = "/home/ajuffer/simploce/pt-cgmd/particles/resources/particles-specs.dat";
    std::ifstream stream;
    file::open_input(stream, fileName);
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
    std::cout << "simulation-model-factory-test test 3" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% model-factory-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    //std::cout << "%TEST_STARTED% test1 (model-factory-test)" << std::endl;
    //test1();
    //std::cout << "%TEST_FINISHED% time=0 test1 (model-factory-test)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (model-factory-test)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (model-factory-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

