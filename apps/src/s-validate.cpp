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
#include "simploce/particle/particle-system.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/util/specification.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/program-options.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>

using namespace simploce;
namespace po = boost::program_options;

static void validate(const p_system_ptr_t& particleSystem,
                     const spec_catalog_ptr_t& catalog,
                     const param_ptr_t& param,
                     const bc_ptr_t& bc) {
    std::cout << std::endl;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);

    auto SBP = catalog->staticBP();
    std::cout << "Validating particle system:" << std::endl;
    std::cout << "Number of particles: " << std::to_string(particleSystem->numberOfParticles()) << std::endl;
    std::cout << "Number of free particles: " << std::to_string(particleSystem->numberOfFreeParticles()) << std::endl;
    std::cout << "Number of frozen particles: " << std::to_string(particleSystem->numberOfFrozenParticles()) << std::endl;
    std::cout << "Number of particle groups: " << std::to_string(particleSystem->numberOfParticleGroups()) << std::endl;
    auto numberSBP = particleSystem->numberOfSpecifications(SBP);
    std::cout << "Number of boundary particles: " << numberSBP << std::endl;
    std::cout << "Box dimension(s)" << *(particleSystem->box()) << std::endl;
    std::cout << "Protonatable?: " << particleSystem->isProtonatable() << std::endl;
    auto isMesoscale = param->get<bool>("simulation.mesoscale");
    std::cout << "Mesoscale?: " << std::to_string(isMesoscale) << std::endl;
    std::cout << "Total charge: " << particleSystem->charge() << std::endl;
    auto M = properties::dipoleMoment(particleSystem, bc);
    std::cout << "Total dipole moment: " <<  M  << std::endl;
    std::cout << "Norm total dipole moment: " << norm<real_t>(M) << std::endl;
    particleSystem->doWithAll<void>([param, isMesoscale] (const std::vector<p_ptr_t>& all) {
        auto p = properties::linearMomentum(all);
        std::cout << "TOTAL linear momentum vector components: " << p << std::endl;
        std::cout << "TOTAL linear momentum vector length: " << norm<real_t>(p) << std::endl;
        energy_t ekin{0.0};
        mass_t total;
        position_t cm{0.0, 0.0, 0.0};

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
    std::cout << std::endl;
    std::cout << "Validating force field." << std::endl;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << *forceField << std::endl;
}

static void validate(const param_ptr_t& param) {
    std::cout << std::endl;
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
        std::string spec = PARTICLE.spec();
        std::string selectFrom = "Input specification. Select one from '" +
                PARTICLE.spec() + "' (default), '" +
                SIM_PARAMS.spec() + "', '" +
                FORCE_FIELD.spec() + "'.";
        bool pbc_1{false};                               // Apply 1D-PBC.
        std::string direction{"z"};                      // Direction.

        po::options_description usage("Usage");
        util::addStandardOptions(usage);
        usage.add_options() (
            "input-spec",
            po::value<std::string>(&spec),
            selectFrom.c_str()
        )(
            "pbc-1",
            po::value<std::string>(&direction),
            "Apply periodicity in the given direction only. Select one of 'x', y', and 'z'."
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << "Validate model:" << std::endl;
            std::cout << usage << "\n";
            return 0;
        }
        if (vm.count("pbc-1")) {
            pbc_1 = true;
            direction = vm["pbc-1"].as<std::string>();
        }
        util::verbose(vm);

        util::Specification specification(spec);
        logger.info("Input model specification: " + specification.spec());
        if (specification == PARTICLE ) {
            if (util::isMesoscale(vm)) {
                param->put("simulation.mesoscale", true);
            } else {
                param->put("simulation.mesoscale", false);
            }
            auto catalog = util::getCatalog(vm);
            auto particleSystem = util::getParticleSystem(vm, param);
            bc_ptr_t bc;
            if (pbc_1) {
                boost::trim(direction);
                bc = factory::pbc1dBB(particleSystem->box(), Direction::valueOf(direction[0]));
            } else {
                 bc = factory::pbc(particleSystem->box());
            }
            validate(particleSystem, catalog, param, bc);
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