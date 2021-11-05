/*
 * Creates particle models.
 * Author: André H. Juffer.
 * Created on 02/11/2021, 10:41.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/protonatable-particle-model-factory.hpp"
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

        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnParticleModel{"particle.system"};
        util::Specification specification = DIATOMIC;
        std::string selectFrom = "Model selection. Choose one from '" +
                DIATOMIC.spec() + "' (molecular oxygen), '" +
                ARGON.spec() + "' (liquid), '" +
                POLARIZABLE_WATER.spec() + "' (coarse grained).";

        po::options_description usage("Usage");
        usage.add_options() (
                "fn-particle-spec-catalog",
                po::value<std::string>(&fnParticleSpecCatalog),
                "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
        )(
                "particle-model-spec",
                po::value<util::Specification>(&specification),
                selectFrom.c_str()
        )(
                "fn-particle-model",
                po::value<std::string>(&fnParticleModel),
                "Output file name for newly created particle model. Default is 'particle.model'."
        )(
                "help", "This help message"
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);
        if (vm.count("help") || argc == 1) {
            std::cout << "Create a particle system:" << std::endl;
            std::cout << usage << "\n";
            return 0;
        }

        logger.info("Selected particle system: " + specification.spec());

        auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
        auto factory = factory::protonatableParticleModelFactory(catalog);
        std::ofstream stream;
        util::open_output_file(stream, fnParticleModel);
        if ( specification == DIATOMIC ) {
            auto diatomic = factory->diatomic(distance_t{0.12}, catalog->O());
            stream << diatomic << std::endl;
        } else if ( specification == POLARIZABLE_WATER ) {
            auto box = factory::box(length_t{7.27});
            auto pWater = factory->polarizableWater(box);
            stream << pWater << std::endl;
        } else if ( specification == ARGON) {
            auto box = factory::box(length_t{3.47786});
            auto argon = factory->argon(box, density_t{1374.0}, temperature_t{94.4});
            stream << argon << std::endl;
        } else {
            util::logAndThrow(logger, specification.spec() + ": No such particle model available.");
        }
        stream.close();
        logger.info("Created particle model stored in '" + fnParticleModel + "'");
    } catch (std::exception& exception) {
        logger.error(exception.what());
    }
}



