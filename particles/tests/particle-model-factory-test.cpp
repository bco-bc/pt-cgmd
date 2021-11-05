/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on October 3, 2019, 1:31 PM
 */

#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/p-types.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void argon(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating liquid argon:" << std::endl;
    particle_model_fact_ptr_t factory = factory::particleModelFactory(catalog);
    box_ptr_t box = factory::box(length_t{3.47786});
    auto atomistic = factory->argon(box);
    std::cout << atomistic << std::endl;
    std::cout << std::endl;
}

void diatomic(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating molecular oxygen:" << std::endl;
    auto spec = catalog->O();
    auto factory = factory::particleModelFactory(catalog);
    Atomistic diatomic = factory->diatomic(0.12, spec);
    std::cout << diatomic << std::endl;
    std::cout << std::endl;
}

void coarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating coarse grained polarizable water model:" << std::endl;
    particle_model_fact_ptr_t factory = factory::particleModelFactory(catalog);
    box_ptr_t box = factory::box(length_t{7.27});
    CoarseGrained coarseGrained = factory->polarizableWater(box);
    std::cout << coarseGrained << std::endl;
    std::cout << std::endl;
}

void electrolyteSolution(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating simple electrolyte solution: " << std::endl;
    auto factory = factory::particleModelFactory(catalog);
    box_ptr_t box = factory::box(length_t{7.0});
    Atomistic electrolyte = factory->simpleElectrolyte(box);
    std::cout << electrolyte << std::endl;
    std::cout << "Number of Na+: ";
    std::cout << boost::lexical_cast<std::string>(electrolyte.numberOfSpecifications(catalog->lookup("Na+")));
    std::cout << std::endl;
}

int main() {
    util::Logger logger("main");
    logger.changeLogLevel(util::Logger::LOGDEBUG);
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fileName);
    std::cout << *catalog << std::endl << std::endl;

    //diatomic(catalog);
    argon(catalog);
    //electrolyteSolution(catalog);
    //coarseGrainedPolarizableWater(catalog);

    return (EXIT_SUCCESS);
}

