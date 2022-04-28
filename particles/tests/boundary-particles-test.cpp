/*
 * Author: André H. Juffer
 *
 * Created on 4/20/22.
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/file.hpp"
#include <iostream>
#include <cstdlib>

using namespace simploce;

spec_catalog_ptr_t getCatalog() {
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    auto catalog = ParticleSpecCatalog::obtainFrom(stream);
    return catalog;
}

p_system_ptr_t getElectrolyte(const spec_catalog_ptr_t& catalog) {
    auto factory = factory::particleSystemFactory(catalog);
    return factory->simpleElectrolyte();
}

int main() {
    auto catalog = getCatalog();
    std::cout << *catalog << std::endl;
    std::cout << std::endl;
    auto electrolyte = getElectrolyte(catalog);
    electrolyte->write(std::cout);
    std::cout << std::endl;

    auto factory = factory::particleSystemFactory(catalog);
    factory->addParticleBoundary(electrolyte, 0.3, Plane::XY);
    electrolyte->write(std::cout);
    std::cout << std::endl;

    return (EXIT_SUCCESS);
}
