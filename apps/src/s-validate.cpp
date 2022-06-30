/*
 * Validates particle systems, etc.
 * Author: André H. Juffer.
 * Created on 03/11/2021, 13:49.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/util/specification.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/program-options.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>

using namespace simploce;
namespace po = boost::program_options;

static void validate(const p_system_ptr_t& particleSystem, const param_ptr_t& param) {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Validating particle system:" << std::endl;
    std::cout << "Number of particles: " << util::toString(particleSystem->numberOfParticles()) << std::endl;
    std::cout << "Box dimension(s)" << *(particleSystem->box()) << std::endl;
    std::cout << "Protonatable?: " << particleSystem->isProtonatable() << std::endl;
    particleSystem->doWithAll<void>([param] (const std::vector<p_ptr_t>& all) {
        auto p = properties::linearMomentum(all);
        std::cout << "TOTAL linear momentum vector components: " << p << std::endl;
        std::cout << "TOTAL linear momentum vector length: " << norm<real_t>(p) << std::endl;
        energy_t ekin{0.0};
        mass_t total;
        position_t cm{0.0, 0.0, 0.0};
        auto isMesoscale = param->get<bool>("simulation.mesoscale");
        for (const auto& particle: all) {
            auto v = particle->velocity();
            auto mass = particle->mass();
            ekin += 0.5 * mass() * inner<real_t>(v, v);
            auto r = particle->position();
            cm += mass() * r;
            total += mass;
        }
        auto temperature = properties::kineticTemperature(all, ekin, isMesoscale);
        std::cout << "Temperature: " << temperature << std::endl;
        cm /= total();
        std::cout << "Center of mass: " << cm << std::endl;
    });
    auto volume = particleSystem->box()->volume();
    auto nd = real_t(particleSystem->numberOfParticles()) / volume;
    std::cout << "Number density: " << nd << std::endl;
}

static void validate(const ff_ptr_t& forceField) {
    std::cout << "Validating force field." << std::endl;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << *forceField << std::endl;
}

static void validate(const param_ptr_t& param) {
    std::cout << "Validating (simulation) parameters." << std::endl;
    param::write(std::cout, *param);
}

int main(int argc, char *argv[]) {
    util::Logger logger{"s-validate-model::main"};

    try {
        // Some default parameters
        param_ptr_t param = factory::simulationParameters();

        const util::Specification PARTICLE{"particle"};
        const util::Specification SIM_PARAMS {"parameters"};
        const util::Specification FORCE_FIELD{"force-field"};

        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnIn{"particle.system"};
        util::Specification specification = PARTICLE;
        std::string selectFrom = "Input specification. Select one from '" +
                PARTICLE.spec() + "' (default), '" +
                SIM_PARAMS.spec() + "', '" +
                FORCE_FIELD.spec() + "'.";

        po::options_description usage("Usage");
        util::addStandardOptions(usage);
        usage.add_options() (
                "input-spec",
                po::value<util::Specification>(&specification),
                selectFrom.c_str()
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << "Validate model:" << std::endl;
            std::cout << usage << "\n";
            return 0;
        }
        util::verbose(vm);

        logger.info("Input model specification: " + specification.spec());
        if (specification == PARTICLE ) {
            auto particleSystem = util::getParticleSystem(vm, param);
            validate(particleSystem, param);
        } else if ( specification == SIM_PARAMS ) {
            util::getParameters(vm, param);
            validate(param);
        } else if (specification == FORCE_FIELD) {
            auto forceField = util::getForceField(vm);
            validate(forceField);
        } else {
            throw std::domain_error(specification.spec() + ": Cannot validate. Not such specification.");
        }

    } catch (std::exception& exception) {
        logger.error(exception.what());
    }
}