/*
 * Validates particle systems, etc.
 * Author: André H. Juffer.
 * Created on 03/11/2021, 13:49.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/util/specification.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/param.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>

using namespace simploce;
namespace po = boost::program_options;

static void validate(const prot_cg_sys_ptr_t& cg) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Validating protonatable coarse grained particle model:" << std::endl;
    std::cout << "Number of particles: " << util::toString(cg->numberOfParticles()) << std::endl;
    std::cout << "Number of beads: " << util::toString(cg->numberOfBeads()) << std::endl;
    std::cout << "Number of protonatable beads: " << util::toString(cg->numberProtonatableBeads()) << std::endl;
    std::cout << "Box dimension(s)" << *(cg->box()) << std::endl;
    std::cout << "Protonatable?: " << cg->isProtonatable() << std::endl;
    cg->doWithAll<void>([] (const std::vector<p_ptr_t>& all) {
        auto p = properties::linearMomentum(all);
        std::cout << "TOTAL linear momentum: " << p << std::endl;
    });
    std::cout << std::endl;
}

static void validate(const at_sys_ptr_t& atomistic) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Validating atomistic particle model:" << std::endl;
    std::cout << "Number of particles: " << util::toString(atomistic->numberOfParticles()) << std::endl;
    std::cout << "Number of atoms: " << util::toString(atomistic->numberOfAtoms()) << std::endl;
    std::cout << "Box dimension(s)" << *(atomistic->box()) << std::endl;
    std::cout << "Protonatable?: " << atomistic->isProtonatable() << std::endl;
    atomistic->doWithAll<void>([] (const std::vector<p_ptr_t>& all) {
        auto p = properties::linearMomentum(all);
        std::cout << "TOTAL linear momentum: " << p << std::endl;
    });
}

static void validate(const ff_ptr_t& forceField) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Validating force field:" << std::endl;
    std::cout << *forceField << std::endl;
}

int main(int argc, char *argv[]) {
    util::Logger logger{"s-validate-model::main"};

    try {
        const util::Specification PARTICLE{"particle"};
        const util::Specification SIM_PARAMS {"simulation-parameters"};
        const util::Specification FORCE_FIELD{"force-field"};

        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnIn{"particle.system"};
        util::Specification specification = PARTICLE;
        std::string selectFrom = "Input specification. Select one from '" +
                PARTICLE.spec() + "' (default), '" + SIM_PARAMS.spec() +
                "', '" + FORCE_FIELD.spec() + "'.";
        bool isCoarseGrained{false};

        po::options_description usage("Usage");
        usage.add_options() (
                "fn-in,i",
                po::value<std::string>(&fnIn),
                "Input file name."
        )(
                "fn-particle-spec-catalog,s",
                po::value<std::string>(&fnParticleSpecCatalog),
                "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
        )(
                "input-spec",
                po::value<util::Specification>(&specification),
                selectFrom.c_str()
        )(
                "coarse-grained,c",
                "Input is a coarse-grained description."
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

        if (specification == PARTICLE ) {
            auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
            std::ifstream stream;
            util::open_input_file(stream, fnIn);
            if ( isCoarseGrained ) {
                logger.debug("Assuming (protonatable) coarse grained particle model.");
                auto particleModel = prot_cg_sys_t::obtainFrom(stream, catalog);
                validate(particleModel);
            } else {
                logger.debug("Assuming atomistic particle model.");
                auto atomistic = Atomistic::obtainFrom(stream, catalog);
                validate(atomistic);
            }
        } else if ( specification == SIM_PARAMS ) {
            std::ifstream stream;
            util::open_input_file(stream, fnIn);
            param_t param;
            param::read(stream, param);
            stream.close();
            param::write(std::cout, param);
        } else if (specification == FORCE_FIELD) {
            auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
            auto forceField = factory::forceField(fnIn, catalog);
            validate(forceField);
        } else {
            throw std::domain_error(specification.spec() + ": Cannot validate.");
        }

    } catch (std::exception& exception) {
        logger.error(exception.what());
    }
}