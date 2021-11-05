/*
 * Author: André H. Juffer.
 * Created on 03/11/2021, 13:49.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-types.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/util/specification.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/util.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>

using namespace simploce;
namespace po = boost::program_options;

static void validate(const prot_cg_mod_ptr_t& cg) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Validating protonatable coarse grained particle model:" << std::endl;
    std::cout << "Number of particles: " << util::toString(cg->numberOfParticles()) << std::endl;
    std::cout << "Number of beads: " << util::toString(cg->numberOfBeads()) << std::endl;
    std::cout << "Number of protonatable beads: " << util::toString(cg->numberProtonatableBeads()) << std::endl;
    std::cout << "Box dimension(s)" << *(cg->box()) << std::endl;
    std::cout << "Protonatable?: " << cg->isProtonatable() << std::endl;
    std::cout << std::endl;
}

static void validate(const at_mod_ptr_t& atomistic) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Validating atomistic particle model:" << std::endl;
    std::cout << "Number of particles: " << util::toString(atomistic->numberOfParticles()) << std::endl;
    std::cout << "Number of atoms: " << util::toString(atomistic->numberOfAtoms()) << std::endl;
    std::cout << "Box dimension(s)" << *(atomistic->box()) << std::endl;
    std::cout << "Protonatable?: " << atomistic->isProtonatable() << std::endl;
}

int main(int argc, char *argv[]) {
    util::Logger logger{"s-validate-model::main"};

    try {
        const util::Specification PARTICLE{"particle"};

        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnIn{"particle.system"};
        util::Specification specification = PARTICLE;
        std::string selectFrom = "Model type. Select one from '" +
                PARTICLE.spec() + "' ...";
        bool isCoarseGrained{false};

        po::options_description usage("Usage");
        usage.add_options() (
                "fn-in",
                po::value<std::string>(&fnIn),
                "Input file name."
        )(
                "fn-particle-spec-catalog",
                po::value<std::string>(&fnParticleSpecCatalog),
                "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
        )(
                "model-specification",
                po::value<util::Specification>(&specification),
                selectFrom.c_str()
        )(
                "coarse-grained",
                "Input is a coarse-grained description."
        )(
                "verbose",
                "Verbose"
        )(
                "help", "This help message"
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);
        if (vm.count("help") || argc == 1) {
            std::cout << "Validate model:" << std::endl;
            std::cout << usage << "\n";
            return 0;
        }
        if (vm.count("coarse-grained") ) {
            isCoarseGrained = true;
        }
        if (vm.count("verbose") ) {
            logger.changeLogLevel(util::Logger::LOGTRACE);
        }

        logger.info("Input model specification: " + specification.spec());

        auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
        auto factory = factory::protonatableParticleModelFactory(catalog);

        if (specification == PARTICLE ) {
            std::ifstream stream;
            util::open_input_file(stream, fnIn);
            if ( isCoarseGrained ) {
                logger.debug("Assuming (protonatable) coarse grained particle model.");
                auto particleModel = prot_cg_mod::obtainFrom(stream, catalog);
                validate(particleModel);
            } else {
                logger.debug("Assuming atomistic particle model.");
                auto atomistic = Atomistic::obtainFrom(stream, catalog);
                validate(atomistic);
            }
        }

    } catch (std::exception& exception) {
        logger.error(exception.what());
    }
}