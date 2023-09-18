/*
 * Author: André H. Juffer.
 * Created on 7/14/2023
 *
 * Copyright (c) 2023 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/s-util.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include <iostream>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGINFO);
    util::Logger logger = util::Logger("one-water-dpd-particle-simulation::main");

    std::ofstream trajectory, data;
    util::open_output_file(trajectory, "/wrk3/simulation/one_water_dpd/trajectory.dat");
    util::open_output_file(data, "/wrk3/simulation/one_water_dpd/simulation.dat");
    logger.info("/wrk3/simulation/one_water_dpd/trajectory.dat: trajectory.");
    logger.info("/wrk3/simulation/one_water_dpd/simulation.dat: simulation data.");

    auto param = factory::simulationParameters();
    param->put("simulation.nsteps", 30000);
    param->put("simulation.npairlists", 10);
    param->put("simulation.nwrite", 20);
    param->put("simulation.timestep", 0.01);
    param->put("simulation.temperature", 1.0);
    param->put("simulation.displacer.dpd.gamma", 2.25);
    param->put("simulation.displacer.dpd.lambda", 0.50);
    param->put("simulation.displacer.dpd.weight-factor", 0.5);
    param->put("simulation.forces.cutoff", 1.0);
    param->put("simulation.mesoscale", 1);
    param->put("simulation.forces.exclude-frozen", true);

    // Thermostat only.
    // param->put("simulation.forces.conservative", false);

    std::clog << "Parameters:" << std::endl;
    param::write(std::clog, param);

    auto catalog = factory::particleSpecCatalog("/wrk3/simulation/one_water_dpd/water-specs.dat");
    auto factory = factory::particleSystemFactory(catalog);
    auto particleSystem = factory::particleSystem("/wrk3/simulation/one_water_dpd/one-water.ps", catalog, true);
    temperature_t temperature = param->get<real_t>("simulation.temperature");
    particleSystem->doWithAll<void>([temperature] (const std::vector<p_ptr_t>& particles) {
        for (auto p: particles) {
            util::assignVelocity(p, temperature, true);
        }
    });
    util::removeOverallLinearMomentum(particleSystem);
    logger.info(std::to_string(particleSystem->numberOfParticles()) + ": Number of particles.");

    auto forceField =
        factory::forceField("/wrk3/simulation/one_water_dpd/water-interaction-parameters-DPD.dat", catalog);
    catalog->write(std::clog);
    std::clog << std::endl;
    logger.info("/wrk3/simulation/one_water_dpd/water-interaction-parameters-DPD.dat: force field parameters");

    //auto bc = factory::no_bc();
    auto bc = factory::pbc(particleSystem->box());
    auto interactor = factory::interactor(param, forceField, bc);
    auto dpdUnits = factory::dpdUnits(1.0, 1.0, 1.0);
    auto displacer = factory::displacer(conf::DPD, param, interactor, bc, dpdUnits);

    Simulation simulation{param, particleSystem, displacer, bc};
    simulation.perform(trajectory, data);
    trajectory.flush();
    data.flush();
    trajectory.close();
    data.close();

    std::ofstream ostream;
    std::string fn = "/wrk3/simulation/one_water_dpd/one-water-out.ps";
    util::open_output_file(ostream, fn);
    ostream << *particleSystem << std::endl;
    ostream.flush();
    ostream.close();

    return (EXIT_SUCCESS);
}