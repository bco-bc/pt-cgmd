/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on October 3, 2019, 1:31 PM
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/p-types.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/file.hpp"
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

p_system_ptr_t coarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating coarse grained polarizable water model:" << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    auto coarseGrained = factory->polarizableWater();
    //std::cout << *coarseGrained << std::endl;
    std::cout << std::endl;
    return coarseGrained;
}

void electrolyteSolution(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating simple electrolyte solution: " << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    box_ptr_t box = factory::box(length_t{7.0});
    auto electrolyte = factory->simpleElectrolyte(box);
    std::cout << *electrolyte << std::endl;
    std::cout << "Number of Na+: ";
    std::cout << util::to_string(electrolyte->numberOfSpecifications(catalog->lookup("Na+")));
    std::cout << std::endl;
}

p_system_ptr_t identicalParticles(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating identical particles:" << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    box_ptr_t box = factory::box(5, 5, 10);
    std::string specName{"H2Om"};
    number_density_t rho{3.0};
    auto particleSystem = factory->identicalParticles(box, specName, rho, temperature_t{1.0});
    std::cout << "Number of particles: " << particleSystem->numberOfParticles() << std::endl;
    //ParticleSystem::validate(particleSystem);
    std::cout << std::endl;
    return particleSystem;
}

p_system_ptr_t mesoscalePolarizableWater(const spec_catalog_ptr_t& catalog) {
    std::cout << "Creating mesoscale polarizable water:" << std::endl;
    auto factory = factory::particleSystemFactory(catalog);
    auto particleSystem = factory->mesoscalePolarizableWater();
    std::cout << "Number of beads: " << particleSystem->numberOfParticles() << std::endl;
    std::cout << "Number of groups: " << particleSystem->numberOfParticleGroups() << std::endl;
    std::cout << "Number of CW beads: " << particleSystem->numberOfSpecifications(catalog->lookup("CW")) << std::endl;
    std::cout << "Number of DP beads: " << particleSystem->numberOfSpecifications(catalog->lookup("DP")) << std::endl;
    return particleSystem;
}

int main() {
    util::Logger logger("main");
    simploce::util::Logger::changeLogLevel(util::Logger::LOGINFO);
    //std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::string fileName = "/wrk3/simulation/meso-water/water-specs.dat";
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fileName);
    std::cout << *catalog << std::endl << std::endl;

    //diatomic(catalog);
    //argon(catalog);
    //electrolyteSolution(catalog);
    // auto particleSystem = coarseGrainedPolarizableWater(catalog);
    // auto particleSystem = identicalParticles(catalog);
    auto particleSystem = mesoscalePolarizableWater(catalog);

    fileName = "/wrk3/tests/particles.ps";
    std::ofstream stream;
    util::open_output_file(stream, fileName);
    stream << *particleSystem << std::endl;
    stream.close();
    logger.info(fileName + ": Particle system written to this file.");

    return (EXIT_SUCCESS);
}

