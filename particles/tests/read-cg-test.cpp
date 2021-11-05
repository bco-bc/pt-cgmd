/*
 * File:   read-cg-argon.cpp
 * Author: juffer
 *
 * Created on 21 September 2019, 14:53
 */

#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "Particle specifications:" << std::endl;
    std::ifstream stream;
    util::open_input_file(stream,
                    "/localdisk/resources/particles-specs.dat");
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::obtainFrom(stream);
    stream.close();
    std::cout << *catalog << std::endl;
    
    std::cout << "Particle model:" << std::endl;
    util::open_input_file(stream,
                    "/localdisk/resources/coarse-grained-system.dat");
    cg_mod_ptr_t cg = CoarseGrained::obtainFrom(stream, catalog);
    stream.close();
    std::cout << *cg << std::endl;
}


int main(int argc, char** argv) {
    test1();
    return (EXIT_SUCCESS);
}

