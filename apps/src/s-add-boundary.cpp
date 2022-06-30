/*
 * Adds boundary particles to existing particle system.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on April 20, 2022, 4.57 pm.
 */

#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>

namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[]) {
    util::Logger logger{"simploce::s-add-boundary::main"};
    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnInputParticleSystem{"in.ps"};
    std::string fnOutputParticleSystem{"out.ps"};
    bool isCoarseGrained{false};

    po::options_description usage("Usage");
    usage.add_options() (
        "fn-particle-spec-catalog,s",
        po::value<std::string>(&fnParticleSpecCatalog),
        "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
    )(
        "fn-input-particle-system,i",
        po::value<std::string>(&fnInputParticleSystem),
        "Input file name particle system. Default is 'in.ps'."
    )(
        "coarse-grained,c",
        "Input is a coarse-grained description."
    )(
        "fn-output-pdb,o",
        po::value<std::string>(&fnOutputParticleSystem),
        "Output file name PDB."
    )(
        "verbose,v",
        "Verbose"
    )(
         "help,h",
         "Help message"
    );

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, usage), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << "Add boundary particles" << std::endl;
        std::cout << usage << "\n";
        return 0;
    }
    if (vm.count("fn-particle-spec-catalog")) {
        fnParticleSpecCatalog = vm["fn-particle-spec-catalog"].as<std::string>();
    }
    if (vm.count("fn-input-particle-system")) {
        fnInputParticleSystem = vm["fn-input-particle-system"].as<std::string>();
    }
    if (vm.count("coarse-grained") ) {
        isCoarseGrained = true;
    }
    if (vm.count("fn-output-particle-system")) {
        fnOutputParticleSystem = vm["fn-output-particle-system"].as<std::string>();
    }
    if (vm.count("verbose") ) {
        util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    }

    // Read particle specifications.
    auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);

    // Read particle system.
    auto particleSystem = factory::particleSystem(fnInputParticleSystem, catalog, isCoarseGrained);
    logger.info("Read particle system from '" + fnInputParticleSystem + "'.");

    // Add boundary particles.
    auto factory = simploce::factory::particleSystemFactory(catalog);
    factory->addParticleBoundary(particleSystem, 0.6, Plane::YZ);
    factory->addParticleBoundary(particleSystem, 0.6, Plane::ZX, true);

    // Write particle system.
    std::ofstream stream;
    util::open_output_file(stream, fnOutputParticleSystem);
    stream << *particleSystem << std::endl;
    stream.close();
    logger.info("Wrote particle system with boundary particles to: " + fnOutputParticleSystem);

    return EXIT_SUCCESS;
}