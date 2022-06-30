/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 15:03.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/file.hpp"
#include <iostream>

using namespace simploce;


int main() {
    util::Logger::changeLogLevel(util::Logger::LOGWARN);

    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto factory = factory::particleSystemFactory(catalog);
    auto particleSystem = simploce::factory::particleSystem("/wrk3/tests/polymer-solution.ps",
                                                            catalog,
                                                            true);

    auto bc = factory::boundaryCondition(particleSystem->box());
    auto forceField =
            factory::forceField("/localdisk/resources/interaction-parameters-polymer-solution.dat", catalog);
    std::clog << *forceField << std::endl;

    auto param = factory::simulationParameters();
    param->put("comment", "For MC and DPD");
    param->put("simulation.nsteps", 20000);
    param->put("simulation.nwrite", 100);
    param->put("simulation.timestep", 0.01);
    param->put("simulation.temperature", 2.0);
    param->put("simulation.displacer.dpd.gamma", 4.5);
    param->put("simulation.displacer.dpd.lambda", 0.50);

    param->put("simulation.displacer.mc.range", 0.5);

    param->put("simulation.forces.cutoff", 1.0);
    param->put("simulation.mesoscale", 1);

    param->put("simulation.scale-velocities", false);
    param->put("simulation.nscale-velocities", 1);
    param->put("simulation.relative-temperature-difference", 10.0);
    param->put("simulation.remove-com-motion", false);
    param->put("simulation.nremove-com-motion", 1);

    // Thermostat only.
    //param->put("simulation.forces.conservative", false);

    param::write(std::clog, *param);

    auto interactor = factory::interactor(param, forceField, bc);
    units::dpd_ptr_t dpdUnits = factory::dpdUnits(1.0, 1.0, 1.0);
    //auto displacer = factory::displacer(conf::MONTE_CARLO, param, interactor, bc);
    auto displacer = factory::displacer(conf::DPD, param, interactor, bc, dpdUnits);
    //auto displacer = factory::displacer(conf::S1_DPD, param, interactor, bc, dpdUnits);
    //auto displacer = factory::displacer(conf::VELOCITY_VERLET, param, interactor, bc, dpdUnits);
    //auto displacer = factory::displacer(conf::LANGEVIN_VELOCITY_VERLET, param, interactor, bc, dpdUnits);
    Simulation simulation{param, particleSystem, displacer};

    std::ofstream trajectoryStream, dataStream;
    util::open_output_file(trajectoryStream, "/wrk3/tests/trajectory.dat");
    util::open_output_file(dataStream, "/wrk3/tests/simulation.dat");
    simulation.perform(trajectoryStream, dataStream);
    trajectoryStream.flush();
    trajectoryStream.close();
    dataStream.flush();
    dataStream.close();

    std::ofstream ostream;
    util::open_output_file(ostream, "/wrk3/tests/polymer-solution-mvv_dpd.ps");
    ostream << *particleSystem << std::endl;

    return 0;
}