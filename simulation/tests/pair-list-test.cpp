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
 * File:   pair-list-test.cpp
 * Author: ajuffer
 *
 * Created on September 4, 2019, 2:55 PM
 */

#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/stypes.hpp"
#include "simploce/simulation/model-factory.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/simulation/pbc.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "pair-list-test test 1" << std::endl;
    
    using generator_t = ParticlePairListGenerator<Bead>;        
    using particle_ptr_t = generator_t::p_ptr_t;
    using particle_ptr_group_t = generator_t::pg_ptr_t;
        
    std::string fileName = "/home/ajuffer/simploce/particles/resources/particles-specs.dat";
    std::ifstream stream;
    file::open_input(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);
    std::cout << *catalog << std::endl;
    
    sim_model_factory_ptr_t mf = factory::modelFactory(catalog);
    
    box_ptr_t box = std::make_shared<box_t>(5.0);
    cg_sim_model_ptr_t polWater = mf->createPolarizableWater(box);

    std::cout << "Number of polarizable waters: " << polWater->size() << std::endl;
    
    bc_ptr_t bc = factory::pbc(box);
    cg_ppair_list_gen_ptr_t generator = factory::coarseGrainedPairListGenerator(box, bc);
    std::vector<particle_ptr_t> free{};
    std::vector<particle_ptr_t> all{};
    std::vector<particle_ptr_group_t> groups{};
    generator->generate(all, free, groups);
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% pair-list-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (pair-list-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (pair-list-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

