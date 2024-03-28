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
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>


namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[]) {
    util::Logger logger{"simploce::s-add-boundary::main"};
    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnUpdatedParticleSpecCatalog{"updated-particle-spec-catalog.dat"};
    std::string fnInputParticleSystem{"in.ps"};
    std::string fnOutputParticleSystem{"out.ps"};
    bool isCoarseGrained{false};
    bool isMesoscale{false};
    bool channel{false};
    bool surface{false};
    real_t boundaryWidth{0.0};
    real_t spacing{0.3};
    real_t sigma{0.0};
    real_t temperature{298.15};
    bool rough{false};
    bool noBoundaryParticles{false};

    po::options_description usage("Usage");
    usage.add_options() (
        "fn-particle-spec-catalog,s",
        po::value<std::string>(&fnParticleSpecCatalog),
        "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
    )(
        "fn-updated-particle-spec-catalog,u",
        po::value<std::string>(&fnUpdatedParticleSpecCatalog),
        "Output file name of updated particle specifications. Default 'updated-particle-spec-catalog.dat'."
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
        "surface",
        "Creates a surface consisting of boundary particles located parallel to the XY plane at z=0."
    )(
        "no-sb-particle",
        "The surface is created without actual particles. However, ions are added to neutralize the system."
    )(
        "srf-cg-density",
        po::value<real_t>(&sigma),
        "Surface charge density. Default is 0.0."
    )(
        "boundary-width",
        po::value<real_t>(&boundaryWidth),
        "Width of the boundary particle layer. Default is 0.0 so that particles are placed on the box boundaries."
    )(
        "spacing",
        po::value<real_t>(&spacing),
        "Spacing between boundary particles, the latter are newly created particles added to "
        "the particle system. This is option is not used when the channel option is selected. "
        "Default width is 0.3."
    )(
        "temperature,T",
        po::value<real_t>(&temperature),
        "Temperature. Default is 298.15 K."
    )(
        "make-rough,r",
        "Make the boundary \"rough\", that is boundary particles are not placed on the box boundary, "
        "but displaced relative to it. Not used when the channel option is selected. Default is "
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
    if (vm.count("fn-updated-particle-spec-catalog")) {
        fnUpdatedParticleSpecCatalog = vm["fn-updated-particle-spec-catalog"].as<std::string>();
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
    if (vm.count("surface")) {
        surface = true;
    }
    if (vm.count("no-sb-particle")) {
        noBoundaryParticles = true;
    }
    if (vm.count("srf-cg-density")) {
        sigma = vm["srf-cg-density"].as<real_t>();
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
    // Create a temporal particle specification for a boundary particle.
    auto specBP = ParticleSpec::create("SBP", 0.0, 1.0, spacing/2.0, true, "# Surface boundary particle");
    auto factory = simploce::factory::particleSystemFactory(catalog);
    auto bc = simploce::factory::pbc(particleSystem->box());
    util::placeInsideBox(particleSystem, bc);
    if (channel) {
        factory->makeChannel(particleSystem, boundaryWidth, isMesoscale);
    } else if (surface) {
        if (noBoundaryParticles) {
            auto adjustedSigma = factory->addChargedSurface(particleSystem, sigma, temperature);
            logger.info(util::to_string(adjustedSigma) + ": Adjusted surface charge density.");
        } else {
            bool bothSides{false};
            bool excludeCorner{false};
            factory->addBoundaryParticles(particleSystem,
                                          spacing,
                                          specBP,
                                          sigma,
                                          Plane::XY,
                                          bothSides,
                                          excludeCorner,
                                          temperature,
                                          isMesoscale,
                                          rough,
                                          boundaryWidth);
        }
    } else {
        bool bothSides{true};
        bool excludeCorner{true};
        // Z-direction remains free of boundary particles.
        factory->addBoundaryParticles(particleSystem,
                                      spacing,
                                      specBP,
                                      0.0,
                                      Plane::YZ,
                                      bothSides,
                                      excludeCorner,
                                      temperature,
                                      isMesoscale,
                                      rough,
                                      boundaryWidth);
        excludeCorner = true;
        factory->addBoundaryParticles(particleSystem,
                                      spacing,
                                      specBP,
                                      0.0,
                                      Plane::ZX,
                                      bothSides,
                                      excludeCorner,
                                      temperature,
                                      isMesoscale,
                                      rough,
                                      boundaryWidth);
    }

    // Write particle system.
    std::ofstream stream;
    util::open_output_file(stream, fnOutputParticleSystem);
    stream << *particleSystem << std::endl;
    stream.close();
    logger.info(fnOutputParticleSystem + ": Particle system with boundary particles.");

    auto specs = particleSystem->specsInUse();
    auto cat = ParticleSpecCatalog::create(specs);
    util::open_output_file(stream, fnUpdatedParticleSpecCatalog);
    stream << *cat << std::endl;
    stream.close();
    logger.info(fnUpdatedParticleSpecCatalog + ": Created new particle specification catalog for particle system.");

    return EXIT_SUCCESS;
}