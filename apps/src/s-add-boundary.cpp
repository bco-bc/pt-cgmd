/*
 * Adds boundary particles to existing particle system.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on April 20, 2022, 4.57 pm.
 */

#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-util.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/box.hpp"
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
    bool isMesoscale{false};
    bool channel{false};
    real_t boundaryWidth{1.0};
    real_t spacing{1.0};
    real_t temperature{298.15};
    bool rough{false};

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
        "fn-output-particle-system,o",
        po::value<std::string>(&fnOutputParticleSystem),
        "Output file name particle system. Default is 'out.ps'."
    )(
        "is-coarse-grained,c",
        "Input is a coarse-grained particle system."
    )(
        "is-mesoscale,m",
        "Input is a mesoscopic particle system."
    )(
        "channel",
        "Creates a channel in the z-direction by creating a wall of particles of given width. "
        "Boundary particles are taken from the given particle system. No additional particles are created."
    )(
        "boundary-width",
        po::value<real_t>(&boundaryWidth),
        "Width of channel boundary. Default is 1.0."
    )(
        "spacing",
        po::value<real_t>(&spacing),
        "Spacing between boundary particles, the latter are newly created particles added to "
        "the particle system. This is option is not applicable when the channel option is selected. "
        "Default width is 1.0."
    )(
        "temperature,T",
        po::value<real_t>(&temperature),
        "Temperature. Default is 298.15 K."
    )(
        "make-rough,r",
        "Make the boundary \"rough\", that is boundary particles are not placed on the box boundary, "
        "but displaced relative to it. Not applicable when the channel option is selected. Default is "
        "to place boundary particle on the box boundaries."
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
    if (vm.count("fn-output-particle-system")) {
        fnOutputParticleSystem = vm["fn-output-particle-system"].as<std::string>();
    }
    if (vm.count("is-coarse-grained") ) {
        isCoarseGrained = true;
    }
    if (vm.count("is-mesoscale") ) {
        isMesoscale = true;
        isCoarseGrained = true;
    }
    if (vm.count("channel")) {
        channel = true;
    }
    if (vm.count("boundary-width")) {
        boundaryWidth = vm["boundary-width"].as<real_t>();
    }
    if (vm.count("spacing")) {
        spacing = vm["spacing"].as<real_t>();
    }
    if (vm.count("temperature")) {
        temperature = vm["temperature"].as<real_t>();
    }
    if (vm.count("make-rough")) {
        rough = true;
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
    auto bc = simploce::factory::pbc(particleSystem->box());
    util::placeInsideBox(particleSystem, bc);
    if (channel) {
        factory->makeChannel(particleSystem, boundaryWidth, isMesoscale);
    } else {
        // Z-direction remains free of boundary particles.
        factory->addParticleBoundary(particleSystem,
                                     spacing,
                                     Plane::YZ,
                                     false,
                                     temperature,
                                     isMesoscale,
                                     rough,
                                     boundaryWidth);
        factory->addParticleBoundary(particleSystem,
                                     spacing,
                                     Plane::ZX,
                                     rough,
                                     temperature,
                                     isMesoscale,
                                     false,
                                     boundaryWidth);
    }

    // Write particle system.
    std::ofstream stream;
    util::open_output_file(stream, fnOutputParticleSystem);
    stream << *particleSystem << std::endl;
    stream.close();
    logger.info("Wrote particle system with boundary particles to: " + fnOutputParticleSystem);

    return EXIT_SUCCESS;
}