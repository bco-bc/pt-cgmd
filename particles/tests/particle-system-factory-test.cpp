/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on October 3, 2019, 1:31 PM
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/p-types.hpp"
#include "simploce/util/util.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void argon(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating liquid argon:" << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    // auto box = factory::box(length_t{3.47786});
    auto argon = factory->argon();
    std::cout << *argon << std::endl;
    std::cout << std::endl;
}

void diatomic(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating molecular oxygen:" << std::endl;
    auto spec = catalog->O();
    auto factory = factory::particleSystemFactory(catalog);
    auto diatomic = factory->diatomic(0.12, spec);
    std::cout << *diatomic << std::endl;
    std::cout << std::endl;
}

void coarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating coarse grained polarizable water model:" << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    auto coarseGrained = factory->polarizableWater();
    std::cout << *coarseGrained << std::endl;
    std::cout << std::endl;
}

void electrolyteSolution(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating simple electrolyte solution: " << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    box_ptr_t box = factory::box(length_t{7.0});
    auto electrolyte = factory->simpleElectrolyte(box);
    std::cout << *electrolyte << std::endl;
    std::cout << "Number of Na+: ";
    std::cout << util::toString(electrolyte->numberOfSpecifications(catalog->lookup("Na+")));
    std::cout << std::endl;
}

int main() {
    util::Logger logger("main");
    logger.changeLogLevel(util::Logger::LOGDEBUG);
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fileName);
    std::cout << *catalog << std::endl << std::endl;

    //diatomic(catalog);
    //argon(catalog);
    electrolyteSolution(catalog);
    //coarseGrainedPolarizableWater(catalog);

    return (EXIT_SUCCESS);
}

