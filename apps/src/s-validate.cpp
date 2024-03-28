/*
 * Validates particle systems, etc.
 * Author: André H. Juffer.
 * Created on 03/11/2021, 13:49.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/util/specification.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/program-options.hpp"
#include <boost/program_options.hpp>
#include <memory>
#include <string>
#include <iostream>

using namespace simploce;
namespace po = boost::program_options;

static void validate(const p_system_ptr_t& particleSystem,
                     const spec_catalog_ptr_t& catalog,
                     bool isMesoscale,
                     const bc_ptr_t& bc,
                     const std::string& sbpSpecName = "SBP") {
    std::cout << std::endl << std::endl;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);

    std::cout << "Validating particle system:" << std::endl;
    std::cout << std::endl;

    spec_ptr_t SBP;
    if (catalog->hasSpecification(sbpSpecName)) {
        SBP = catalog->lookup(sbpSpecName);
    } else {
        SBP = ParticleSpec::create(sbpSpecName, 0.0, 1.0, 1.0, true, "SBP");
    }
    auto box = particleSystem->box();
    auto volume = particleSystem->box()->volume();
    auto specs = particleSystem->specsInUse();
    size_t total = 0;
    for (auto& spec: specs) {
        auto n = particleSystem->numberOfSpecifications(spec.second);
        auto name = spec.second->name();
        std::cout << std::to_string(n) + ": Number of particles of specification '" + name + "'." << std::endl;
        if (name == SBP->name()) {
            // Surface boundary particles.
            std::cout << name  + ": Assuming particles of this specification represent surface boundary particles:" << std::endl;
            auto charge = n * spec.second->charge();
            auto sigma = charge() / (box->lengthX() * box->lengthY());
            std::cout << "   " << std::to_string(sigma) << ": Surface charge density (xy-plane)." << std::endl;
            std::cout << "   " << std::to_string(charge()) << ": Surface total charge." << std::endl;
        } else {
            // Other type of particles.
            auto nd = real_t(n) / volume;
            std::cout << nd << ": Number density '" + name + "'." << std::endl;
            total += n;
        }
    }
    std::cout << "Total charge: " << particleSystem->charge() << std::endl;
    auto M = properties::dipoleMoment(particleSystem, bc);
    std::cout << "Total dipole moment: " <<  M  << std::endl;
    std::cout << "Norm total dipole moment: " << norm<real_t>(M) << std::endl;
    particleSystem->doWithAll<void>([isMesoscale] (const std::vector<p_ptr_t>& all) {
        auto p = properties::linearMomentum(all);
        std::cout << p << ": TOTAL linear momentum vector components: " << std::endl;
        std::cout << norm<real_t>(p) <<": TOTAL linear momentum vector length: " << std::endl;
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
        std::cout << temperature << ": Temperature" << std::endl;
        cm /= total();
        std::cout << cm << ": Center of mass" << std::endl;
    });
    auto nd = real_t(particleSystem->numberOfParticles()) / volume;
    std::cout << nd << ": TOTAL Number density: " << std::endl;
}

static void validate(const ff_ptr_t& forceField) {
    std::cout << std::endl << std::endl;
    std::cout << "Validating force field." << std::endl;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << *forceField << std::endl;
}

static void validate(const param_ptr_t& param) {
    std::cout << std::endl << std::endl;
    std::cout << "Validating (simulation) parameters." << std::endl;
    param::write(std::cout, *param);
}

int main(int argc, char *argv[]) {
    util::Logger logger{"s-validate::main"};

    try {
        const util::Specification PARTICLE{"particle"};
        const util::Specification SIM_PARAMS {"parameters"};
        const util::Specification FORCE_FIELD{"force-field"};

        std::string spec = PARTICLE.spec();
        std::string selectFrom = "Input specification. Select one from '" +
                PARTICLE.spec() + "' (default), '" +
                SIM_PARAMS.spec() + "', '" +
                FORCE_FIELD.spec() + "'.";
        bool pbc_1{false};                               // Whether to apply 1D-PBC.
        bool pbc_2{false};                               // Whether to apply 2D-PBC.
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
        )(
            "pbc-2",
            "Assume periodicity in x- and y-direction, but not in the z-direction."
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
        if (vm.count("pbc-2")) {
            pbc_2 = true;
        }
        util::verbose(vm);

        util::Specification specification(spec);
        logger.info("Input model specification: " + specification.spec());
        if (specification == PARTICLE ) {
            param_ptr_t param = std::make_shared<param_t>();  // Empty.
            auto catalog = util::getCatalog(vm);
            auto particleSystem = util::getParticleSystem(vm, param);
            bc_ptr_t bc;
            if (pbc_1) {
                boost::trim(direction);
                bc = factory::pbc1dBB(particleSystem->box(), Direction::valueOf(direction[0]));
            } else if (pbc_2) {
                bc = factory::pbc_2d(particleSystem->box(), Direction::X, Direction::Y);
            } else {
                 bc = factory::pbc(particleSystem->box());
            }
            bool isMesoscale = param->get<bool>("simulation.mesoscale");
            validate(particleSystem, catalog, isMesoscale, bc);
        } else if ( specification == SIM_PARAMS ) {
            param_ptr_t param = std::make_shared<param_t>();
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