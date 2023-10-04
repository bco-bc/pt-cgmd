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
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/file.hpp"
#include <iostream>

using namespace simploce;


int main() {
    util::Logger::changeLogLevel(util::Logger::LOGWARN);

    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    auto factory = factory::particleSystemFactory(catalog);


    auto particleSystem =
            simploce::factory::particleSystem("/wrk3/tests/particles-in-channel.ps",
                                              catalog,
                                              true);


    /*
    auto particleSystem =
            simploce::factory::particleSystem("/wrk3/tests/identical-particles.ps",
                                              catalog,
                                              true);
    */

    util::Logger::changeLogLevel(util::Logger::LOGWARN);
    std::cout << "Number of particles: " << particleSystem->numberOfParticles() << std::endl;
    auto spec = catalog->staticBP();
    std::cout << "Number of boundary particles: " << particleSystem->numberOfSpecifications(spec) << std::endl;
    std::cout << "Number of frozen particles (before freezing): " << particleSystem->numberOfFrozenParticles() << std::endl;
    particleSystem->freeze(spec);
    std::cout << "Number of frozen particles (after freezing): " << particleSystem->numberOfFrozenParticles() << std::endl;
    //util::Logger::changeLogLevel(util::Logger::LOGWARN);

    //auto bc = factory::pbc(particleSystem->box());
    auto bc = factory::pbc1dBB(particleSystem->box(), Direction::Z);
    auto forceField =
            factory::forceField("/localdisk/resources/interaction-parameters-droplets-polymer-solution.dat",
                                catalog);
    std::clog << *forceField << std::endl;

    auto param = factory::simulationParameters();
    param->put("comment", "For MC and DPD");
    param->put("simulation.nsteps", 10000);
    param->put("simulation.npairlists", 10);
    param->put("simulation.nwrite", 100);
    param->put("simulation.timestep", 0.01);
    param->put("simulation.temperature", 1.0);
    param->put("simulation.displacer.dpd.gamma", 2.25);
    param->put("simulation.displacer.dpd.lambda", 0.50);
    param->put("simulation.displacer.dpd.weight-factor", 0.5);

    param->put("simulation.displacer.mc.range", 0.5);

    param->put("simulation.forces.exclude-frozen", true);
    param->put("simulation.forces.include-external", true);
    param->put("simulation.forces.cutoff", 1.0);
    param->put("simulation.mesoscale", 1);

    param->put("simulation.scale-velocities", false);
    param->put("simulation.nscale-velocities", 100);
    param->put("simulation.relative-temperature-difference", 10.0);
    param->put("simulation.remove-com-motion", false);
    param->put("simulation.nremove-com-motion", 1);

    // Thermostat only.
    param->put("simulation.forces.conservative", false);

    std::clog << "Parameters:" << std::endl;
    param::write(std::clog, param);

    auto interactor = factory::interactor(param, forceField, bc);
    units::dpd_ptr_t dpdUnits = factory::dpdUnits(1.0, 1.0, 1.0);
    //auto displacer = factory::displacer(conf::MONTE_CARLO, param, interactor, bc);
    auto displacer = factory::displacer(conf::DPD, param, interactor, bc, dpdUnits);
    //auto displacer = factory::displacer(conf::S1_DPD, param, interactor, bc, dpdUnits);
    //auto displacer = factory::displacer(conf::VELOCITY_VERLET, param, interactor, bc, dpdUnits);
    //auto displacer = factory::displacer(conf::LANGEVIN_VELOCITY_VERLET, param, interactor, bc, dpdUnits);
    Simulation simulation{param, particleSystem, catalog, displacer, bc};

    //util::Logger::changeLogLevel(util::Logger::LOGTRACE);
    std::ofstream trajectoryStream, dataStream;
    util::open_output_file(trajectoryStream, "/wrk3/tests/trajectory.dat");
    util::open_output_file(dataStream, "/wrk3/tests/simulation.dat");
    //util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    simulation.perform(trajectoryStream, dataStream);
    //std::clog << "Parameters (3):" << std::endl;
    //util::Logger::changeLogLevel(util::Logger::LOGTRACE);
    trajectoryStream.flush();
    trajectoryStream.close();
    dataStream.flush();
    dataStream.close();

    std::ofstream ostream;
    std::string fn = "/wrk3/tests/particles-in-channel-dpd.ps";
    util::open_output_file(ostream, fn);
    ostream << *particleSystem << std::endl;
    ostream.flush();
    ostream.close();

    return 0;
}