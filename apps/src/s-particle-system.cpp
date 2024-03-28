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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace po = boost::program_options;
using namespace simploce;

/**
 * Creates a box.
 * @param boxSizes
 * @return Box
 */
static box_ptr_t
makeBox(const std::string& boxSizes) {
    std::vector<std::string> result;
    boost::split(result, boxSizes, boost::is_any_of(","));
    if (result.size() != 3) {
        throw std::domain_error(boxSizes + ": No box sizes. No three box sizes were specified.");
    }
    real_t sizes[3];
    for (int k = 0; k != 3; ++k) {
        auto size = result[k];
        boost::trim(size);
        sizes[k] = boost::lexical_cast<real_t>(size);
    }
    return factory::box(sizes[0], sizes[1], sizes[2]);
}

int main(int argc, char *argv[]) {
    util::Logger logger{"s-particle-system::main"};

    try {
        int numberOfFluidElements = 6000;                      // Number of water fluid elements.
        real_t characteristicLengthMesoscopicWater = 6.0e-06;  // In m.
        real_t boxSize(6.30);                                  // nm.
        std::string boxSizes{"20.0,20.0,40.0"};                // Box sizes.
        real_t molarity{0.1};                                  // Electrolyte molarity (mol/l)

        const util::Specification DIATOMIC{"diatomic"};
        const util::Specification ARGON{"argon"};
        const util::Specification CG_POLARIZABLE_WATER{"cg-polarizable-water"};
        const util::Specification MESOSCOPIC_POLARIZABLE_WATER{"meso-polarizable-water"};
        const util::Specification ELECTROLYTE{"electrolyte"};
        const util::Specification POLYMER_SOLUTION{"polymer-solution"};
        const util::Specification DROPLETS_POLYMER_SOLUTION{"channel-droplets-polymer-solution"};
        const util::Specification WATER_IN_MICROCHANNEL("water-in-microchannel");
        const util::Specification LARGE_OBJECT_IN_ELECTROLYTE("large-object-in-electrolyte");
        bool addBoundaryParticles{false};

        real_t temperature{units::si<real_t>::ROOM_TEMPERATURE};

        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnUpdatedParticleSpecCatalog{"particle-specs-catalog.dat"};
        std::string fnParticleSystem{"out.ps"};
        std::string spec = DIATOMIC.spec();
        std::string selectFrom =
            "Particle system selection. Choose one from '" +
            DIATOMIC.spec() + "' (molecular oxygen, atomistic), '" +
            ARGON.spec() + "' (atomistic, default), '" +
            ELECTROLYTE.spec() + "' (electrolyte), '" +
            LARGE_OBJECT_IN_ELECTROLYTE.spec() + "' (large object in an electrolyte), '" +
            POLYMER_SOLUTION.spec() + "' (polymer solution), '" +
            CG_POLARIZABLE_WATER.spec() + "' (coarse grained),  '" +
            MESOSCOPIC_POLARIZABLE_WATER.spec() + "' (mesoscopic), '" +
            DROPLETS_POLYMER_SOLUTION.spec() + "' (droplets in a polymer solution in a channel), '" +
                WATER_IN_MICROCHANNEL.spec() + "' (water in a microchannel)'.";

        po::options_description usage("Usage");
        usage.add_options() (
            "fn-particle-spec-catalog,s",
            po::value<std::string>(&fnParticleSpecCatalog),
            "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
        )(
            "particle-system-spec,p",
            po::value<std::string>(&spec),
            selectFrom.c_str()
        )(
            "fn-updated-particle-spec-catalog,u",
            po::value<std::string>(&fnUpdatedParticleSpecCatalog),
            "Output file name of updated particle specifications. Default 'updated-particle-spec-catalog.dat'."
        )(
            "add-boundary-particles,b",
            "Add boundary particles (for \"electrolyte\" only. DEPRECATED: Use 's-add-boundary' instead.)"
        )(
            "fn-particle-system,o",
            po::value<std::string>(&fnParticleSystem),
            "Output file name for newly created particle system. Default is 'out.ps'."
        )(
            "temperature,T",
            po::value<real_t>(&temperature),
            "Temperature. Default is 298.15. For mesoscale, select 1, 2, 3, and so forth."
        )(
            "box-size",
            po::value<real_t>(&boxSize),
            "Boz size in one dimension, to be applied in all dimensions. "
            "Default is 6.3 nm for an electrolyte and a large object in electrolyte."
        )(
            "box-sizes",
            po::value<std::string>(&boxSizes),
            "Comma-separated box sizes for three dimensions."
        )(
            "characteristic-length-water,d",
            po::value<real_t>(&characteristicLengthMesoscopicWater),
            "Characteristic length for \"mesoscopic polarizable water\" model. Default is 6.0e-06 m."
        )(
            "n-water-fluid-elements,M", po::value<int>(&numberOfFluidElements),
            "Number of water fluid elements."
        )(
            "molarity,m",
            po::value<real_t>(&molarity),
            "Molarity in mol/l. For electrolyte only."
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
            logger.warn("Adding boundary particles: DEPRECATED. Use 's-add-boundary' instead.");
        }
        if (vm.count("verbose") ) {
            util::Logger::changeLogLevel(util::Logger::LOGTRACE);
        }

        util::Specification specification(spec);
        logger.info("Selected particle system: " + specification.spec());
        logger.info(std::to_string(temperature) + ": Temperature (arbitrary units).");

        auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
        auto factory = factory::protonatableParticleSystemFactory(catalog);
        std::ofstream stream;
        util::open_output_file(stream, fnParticleSystem);

        p_system_ptr_t particleSystem{};

        if ( specification == DIATOMIC ) {
            auto diatomic = factory->diatomic(dist_t{0.12}, catalog->O());
            stream << *diatomic << std::endl;
            particleSystem = diatomic;
        } else if (specification == CG_POLARIZABLE_WATER ) {
            auto pWater = factory->cgmdPolarizableWater();
            stream << *pWater << std::endl;
            particleSystem = pWater;
        } else if (specification == MESOSCOPIC_POLARIZABLE_WATER) {
            // Fixed box dimensions for x- and y-direction. From Serhatlioglu2020.
            real_t lx = 6.0e-05 / characteristicLengthMesoscopicWater;  // In DPD units.
            real_t ly = 5.0e-05 / characteristicLengthMesoscopicWater;
            density_t rho = 3.0;  // Fixed total bead number density for DPD, rho = N / (lx * ly *lz), N=2M, and M is number
                                  // fluid elements.
            auto N = 2 * numberOfFluidElements;
            auto lz = N / (lx * ly * rho());
            auto box = factory::box(lx, ly, lz);
            auto pWater =
                    factory->mesoscalePolarizableWater(box,
                                                       numberOfFluidElements,
                                                       temperature);
            stream << *pWater << std::endl;
            particleSystem = pWater;
        } else if (specification == ARGON) {
            auto argon = factory->argon();
            stream << *argon << std::endl;
            particleSystem = argon;
        } else if (specification == ELECTROLYTE) {
            box_ptr_t box;
            if (vm.count("molarity")) {
                molarity = vm["molarity"].as<real_t>();
            }
            if (vm.count("box-size")) {
                boxSize = vm["box-size"].as<real_t>();
                box = factory::box(boxSize, boxSize, boxSize);
            } else if (vm.count("box-sizes")) {
                box = makeBox(vm["box-sizes"].as<std::string>());
            } else {
                box = factory::box(boxSize, boxSize, boxSize);
            }
            auto electrolyte = factory->simpleElectrolyte(box, molarity, temperature, true);
            stream << *electrolyte << std::endl;
            particleSystem = electrolyte;
            if (addBoundaryParticles) {
                logger.warn("Adding boundary particles is ignored. Use 's-add-boundary' instead.");
            }
        } else if (specification == LARGE_OBJECT_IN_ELECTROLYTE) {
            if (vm.count("box-size")) {
                boxSize = vm["box-size"].as<real_t>();
            }
            auto box = factory::box(boxSize, boxSize, boxSize);
            auto loInElectrolyte = factory->largeObjectInElectrolyte(box);
            stream << *loInElectrolyte << std::endl;
            particleSystem = loInElectrolyte;
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
            particleSystem = polymerSolution;
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
            particleSystem = dropletsPolymerSolution;
        } else if (specification == WATER_IN_MICROCHANNEL) {
            auto box = factory::box(40, 42, 80);
            auto identical =
                factory->identicalParticles(box,
                                            "H2Om",
                                            number_density_t{3.0},
                                            temperature_t{1.0},
                                            true);
            logger.info("Writing particle system...");
            stream << *identical << std::endl;
            logger.info("Done.");
            particleSystem = identical;
        } else {
            util::logAndThrow(logger, specification.spec() + ": No such particle model available.");
        }

        stream.close();
        logger.info("Created particle model stored in '" + fnParticleSystem + "'");

        auto specs = particleSystem->specsInUse();
        auto cat = ParticleSpecCatalog::create(specs);
        util::open_output_file(stream, fnUpdatedParticleSpecCatalog);
        stream << *cat << std::endl;
        stream.close();
        logger.info(fnUpdatedParticleSpecCatalog +
                    ": Particle specification catalog for newly created particle system.");
    } catch (std::exception& exception) {
        logger.error(exception.what());
    }
}



