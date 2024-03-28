/*
 * Author: André H. Juffer
 *
 * Created on 4/20/22.
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec.hpp"
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
    auto box = factory::box(20.0, 20.0, 40.0);
    return factory->simpleElectrolyte(box);
}

p_system_ptr_t getWaterInChannel(const spec_catalog_ptr_t& catalog) {
    auto factory = simploce::factory::particleSystemFactory(catalog);
    box_ptr_t box = factory::box(10.00000, 8.33333, 166.66667);
    std::string specName{"H2Om"};
    number_density_t rho{3.0};
    return factory->identicalParticles(box, specName, rho, temperature_t{1.0});
}

int main() {
    util::Logger logger("main");
    util::Logger::changeLogLevel(util::Logger::LOGDEBUG);

    auto catalog = getCatalog();
    std::cout << *catalog << std::endl;
    std::cout << std::endl;
    auto factory = factory::particleSystemFactory(catalog);

    // Electrolyte
    auto electrolyte = getElectrolyte(catalog);
    spec_ptr_t spec = ParticleSpec::create("SBP", 0.1, 1.0, 0.1, true, "# Surface boundary particle");
    factory->addBoundaryParticles(electrolyte, 0.3, spec);
    std::string fileName = "/wrk3/tests/electrolyte-bp.ps";
    std::ofstream ostream;
    util::open_output_file(ostream, fileName);
    electrolyte->write(ostream);
    ostream.close();

    auto charge = electrolyte->charge();
    auto box = electrolyte->box();
    std::cout << charge << ": Total charge." << std::endl;

    std::cout << "Particle specifications in use:" << std::endl;
    auto specs = electrolyte->specsInUse();
    catalog = ParticleSpecCatalog::create(specs);
    std::cout << *catalog << std::endl;
    spec = catalog->lookup("SBP");
    auto n = electrolyte->numberOfSpecifications(spec);
    auto sigma = n * spec->charge() / (box->lengthX() * box->lengthY());
    std::cout << sigma << ": Surface charge density." << std::endl;
    std::cout << electrolyte->numberOfSpecifications(spec) << ": Number of boundary particles." << std::endl;

    fileName = "/wrk3/tests/catalog.dat";
    util::open_output_file(ostream, fileName);
    ostream << *catalog << std::endl;
    ostream.close();
    std::cout << "New particle specification catalog written to '" + fileName + "'." << std::endl;

    /*
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
     */

    return (EXIT_SUCCESS);
}
