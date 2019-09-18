/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particle-spec-catalog-test.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:50 PM
 */

#include "simploce/util/file.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/ptypes.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "particle-spec-catalog-test test 1" << std::endl;
    std::string fileName = "/home/ajuffer/simploce/particles/resources/particles-specs.dat";
    std::ifstream stream;
    file::open_input(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::create(stream);
    std::cout << *catalog << std::endl;
    atom_spec_ptr_t spec = catalog->lookup("O");
    std::cout << "spec: " << *spec << std::endl;
}

void test2() {
    std::cout << "particle-spec-catalog-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (particle-spec-catalog-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% particle-spec-catalog-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (particle-spec-catalog-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (particle-spec-catalog-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test2 (particle-spec-catalog-test)\n" << std::endl;
    //test2();
    //std::cout << "%TEST_FINISHED% time=0 test2 (particle-spec-catalog-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

