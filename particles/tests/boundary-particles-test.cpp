/*
 * Author: André H. Juffer
 *
 * Created on 4/20/22.
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
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

p_system_ptr_t getWaterInChannel(const spec_catalog_ptr_t& catalog) {
    auto factory = simploce::factory::particleSystemFactory(catalog);
    box_ptr_t box = factory::box(10.00000, 8.33333, 166.66667);
    std::string specName{"H2Om"};
    number_density_t rho{3.0};
    return factory->identicalParticles(box, specName, rho, temperature_t{1.0});
}

int main() {
    auto catalog = getCatalog();
    std::cout << *catalog << std::endl;
    std::cout << std::endl;
    /*
    auto electrolyte = getElectrolyte(catalog);
    electrolyte->write(std::cout);
    std::cout << std::endl;

    auto factory = factory::particleSystemFactory(catalog);
    factory->addParticleBoundary(electrolyte, 0.3, Plane::XY);
    electrolyte->write(std::cout);
    std::cout << std::endl;
    */


    util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    auto particleSystem = getWaterInChannel(catalog);
    auto factory = simploce::factory::particleSystemFactory(catalog);
    factory->makeChannel(particleSystem, length_t{0.1}, true);
    std::string fileName = "/wrk3/tests/particles-in-channel.ps";
    std::ofstream stream;
    util::open_output_file(stream, fileName);
    stream << *particleSystem << std::endl;
    stream.close();
    std::clog << fileName << ": Particle system written to this file." << std::endl;
    auto spec = catalog->staticBP();
    auto numberOfBoundaryParticles = particleSystem->numberOfSpecifications(spec);
    particleSystem->freeze(spec);
    std::clog << numberOfBoundaryParticles << ": Number of boundary particles." << std::endl;
    std::clog << particleSystem->numberOfParticles() << ": Total number of particles." << std::endl;
    return (EXIT_SUCCESS);
}
