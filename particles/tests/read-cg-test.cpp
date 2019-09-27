/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   read-cg-test.cpp
 * Author: juffer
 *
 * Created on 21 September 2019, 14:53
 */

#include "simploce/particle/pall.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "read-cg-test test 1" << std::endl;
    std::cout << "Particle specifications:" << std::endl;
    std::ifstream stream;
    file::open_input(stream, 
                    "/home/ajuffer/simploce/pt-cgmd/particles/resources/particles-specs.dat");
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);
    stream.close();
    std::cout << *catalog << std::endl;
    
    std::cout << "Particle model:" << std::endl;
    file::open_input(stream,
                    "/home/ajuffer/simploce/pt-cgmd/particles/resources/coarse-grained-system.dat");
    cg_ptr_t cg = CoarseGrained::readFrom(stream, catalog);
    stream.close();
    std::cout << *cg << std::endl;
}


int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% read-cg-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (read-cg-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (read-cg-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

