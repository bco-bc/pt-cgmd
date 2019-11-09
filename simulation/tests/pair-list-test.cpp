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

#include "simploce/simulation/cell-lists.hpp"
#include "simploce/simulation/distance-lists.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/stypes.hpp"
#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/sim-util.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "pair-list-test test 1" << std::endl;
    
    using cell_generator_t = CellLists<Bead>;   
    using dist_generator_t = DistanceLists<Bead>; 
    using pairlists_t = PairLists<Bead>;
        
    std::string fileName = "/home/ajuffer/simploce/pt-cgmd/particles/resources/particles-specs.dat";
    std::ifstream stream;
    file::open_input(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);
    //std::clog << *catalog << std::endl;
    
    auto box = factory::cube(7.27);
    auto smf = factory::simulationModelFactory(catalog);
    auto sm = smf->polarizableWater(box, 997.0479, 298.15, 2560);
    
    std::clog << "Number of beads: " << sm->size() << std::endl;
    std::clog << std::endl;
    
    dist_generator_t distGen(sm->box(), sm->boundaryCondition());

    auto start = std::chrono::system_clock::now();
    auto pairlists = 
        sm->doWithAllFreeGroups<pairlists_t>([distGen] (const std::vector<bead_ptr_t>& all,
                                                        const std::vector<bead_ptr_t>& free,
                                                        const std::vector<bead_group_ptr_t>& groups) {
        return distGen.generate(all, free, groups);
    });
    auto end = std::chrono::system_clock::now();
    auto distDiff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::clog << "Distance-based pair lists: Time used: " << distDiff << std::endl;   
    std::clog << std::endl;

    cell_generator_t cellGen(sm->box(), sm->boundaryCondition());
    //cellGen.sideLength(0.25 * util::cutoffDistance(box));
    start = std::chrono::system_clock::now();
    pairlists = 
        sm->doWithAllFreeGroups<pairlists_t>([cellGen] (const std::vector<bead_ptr_t>& all,
                                                        const std::vector<bead_ptr_t>& free,
                                                        const std::vector<bead_group_ptr_t>& groups) {
        return cellGen.generate(all, free, groups);
    });
    end = std::chrono::system_clock::now();
    auto cellDiff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::clog << "Cell-based pair lists: Time used: " << cellDiff << std::endl;  
    std::clog << "Speed ratio (distance/cell): " << real_t(distDiff)/real_t(cellDiff)
              << std::endl;
    std::clog << std::endl;
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

