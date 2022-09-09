/*
 * File:   particle-spec-catalog-argon.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:50 PM
 */

#include "simploce/util/file.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/p-types.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

void test() {
    std::cout << "particle-spec-catalog-test argon 1" << std::endl;
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::obtainFrom(stream);
    std::cout << *catalog << std::endl;
    std::cout << std::endl;
    std::cout << "Lookup O" << std::endl;
    spec_ptr_t spec = catalog->lookup("O");
    std::cout << *spec << std::endl;
    std::cout << "Lookup O" << std::endl;
    spec = catalog->O();
    std::cout << *spec << std::endl;
    std::cout << "Lookup C" << std::endl;
    spec = catalog->C();
    std::cout << *spec << std::endl;
    std::cout << "Lookup N" << std::endl;
    spec = catalog->N();
    std::cout << *spec << std::endl;
}

int main() {
    test();
    return (EXIT_SUCCESS);
}

