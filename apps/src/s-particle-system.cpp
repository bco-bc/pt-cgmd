/*
 * Creates particle systems.
 * Author: André H. Juffer.
 * Created on 02/11/2021, 10:41.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/specification.hpp"
#include <boost/program_options.hpp>
#include <iostream>

namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[]) {
    util::Logger logger{"s-particle-system::main"};

    try {
        const util::Specification DIATOMIC{"diatomic"};
        const util::Specification ARGON{"argon"};
        const util::Specification POLARIZABLE_WATER{"polarizable-water"};
        const util::Specification ELECTROLYTE{"electrolyte"};
        const util::Specification POLYMER_SOLUTION{"polymer-solution"};
        const util::Specification DROPLETS_POLYMER_SOLUTION{"channel-droplets-polymer-solution"};
        const util::Specification WATER_IN_MICROCHANNEL("water-in-microchannel");
        bool addBoundaryParticles{false};

        real_t temperature{units::si<real_t>::ROOM_TEMPERATURE};

        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnParticleSystem{"out.ps"};
        util::Specification specification = DIATOMIC;
        std::string selectFrom =
                "Particle system selection. Choose one from '" +
                DIATOMIC.spec() + "' (molecular oxygen, atomistic), '" +
                ARGON.spec() + "' (atomistic, default), '" +
                ELECTROLYTE.spec() + "' (electrolyte), '" +
                POLYMER_SOLUTION.spec() + "' (polymer solution), '" +
                POLARIZABLE_WATER.spec() + "' (coarse grained),  '" +
                DROPLETS_POLYMER_SOLUTION.spec() + "' (droplets in a polymer solution in a channel), '" +
                WATER_IN_MICROCHANNEL.spec() + "' (water in a microchannel)'.";

        po::options_description usage("Usage");
        usage.add_options() (
                "fn-particle-spec-catalog,s",
                po::value<std::string>(&fnParticleSpecCatalog),
                "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
        )(
                "particle-system-spec,p",
                po::value<util::Specification>(&specification),
                selectFrom.c_str()
        )(
                "add-boundary-particles,b",
                "Add boundary particles"
        )(
                "fn-particle-system,o",
                po::value<std::string>(&fnParticleSystem),
                "Output file name for newly created particle system. Default is 'out.ps'."
        )(
                "temperature,T", po::value<real_t>(&temperature),
                "Temperature. Default is 298.15. For mesoscale, select 1, 2, 3, and so forth."
        )(
                "verbose,v",
                "Verbose"
        )(
                "help,h", "This help message"
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);
        if (vm.count("help") || argc == 1) {
            std::cout << "Create a particle system. The particle systems listed below are based on specific ";
            std::cout << "publications and are created according to the article's authors specifications." << std::endl;
            std::cout << usage << "\n";
            return 0;
        }
        if (vm.count("add-boundary-particles")) {
            addBoundaryParticles = true;
        }
        if (vm.count("verbose") ) {
            util::Logger::changeLogLevel(util::Logger::LOGTRACE);
        }

        logger.info("Selected particle system: " + specification.spec());
        logger.info(std::to_string(temperature) + ": Temperature (arbitrary units).");

        auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
        auto factory = factory::protonatableParticleSystemFactory(catalog);
        std::ofstream stream;
        util::open_output_file(stream, fnParticleSystem);
        if ( specification == DIATOMIC ) {
            auto diatomic = factory->diatomic(dist_t{0.12}, catalog->O());
            stream << *diatomic << std::endl;
        } else if ( specification == POLARIZABLE_WATER ) {
            auto pWater = factory->polarizableWater();
            stream << *pWater << std::endl;
        } else if ( specification == ARGON) {
            auto argon = factory->argon();
            stream << *argon << std::endl;
        } else if ( specification == ELECTROLYTE) {
            auto electrolyte = factory->simpleElectrolyte();
             if (addBoundaryParticles) {
                 factory->addParticleBoundary(electrolyte, 0.6, Plane::YZ);
                 factory->addParticleBoundary(electrolyte, 0.6, Plane::ZX);
             }
            stream << *electrolyte << std::endl;
        } else if ( specification == POLYMER_SOLUTION ) {
            auto box = factory::box(10, 10, 40);
            auto polymerSolution =
                    factory->polymerSolution(box,
                                             15,
                                             "PMU",
                                             400,
                                             1.0,
                                             6000,
                                             "H2Om",
                                             temperature,
                                             true);
            stream << *polymerSolution << std::endl;
        } else if ( specification == DROPLETS_POLYMER_SOLUTION) {
            auto box = factory::box(15, 15, 30);
            auto dropletsPolymerSolution =
                    factory->dropletPolymerSolution(box,
                                            1,
                                            "PMU",
                                            0,
                                            1.0,
                                            10125,
                                            "H2Om",
                                            10125,
                                            "DROm",
                                            1.0);
            stream << *dropletsPolymerSolution << std::endl;
        } else if (specification == WATER_IN_MICROCHANNEL) {
            auto box = factory::box(40, 42, 80);
            auto particleSystem =
                factory->identicalParticles(box,
                                            "H2Om",
                                            number_density_t{3.0},
                                            temperature_t{1.0},
                                            true);
            logger.info("Writing particle system...");
            stream << *particleSystem << std::endl;
            logger.info("Done.");
        } else {
            util::logAndThrow(logger, specification.spec() + ": No such particle model available.");
        }

        stream.close();
        logger.info("Created particle model stored in '" + fnParticleSystem + "'");
    } catch (std::exception& exception) {
        logger.error(exception.what());
    }
}



