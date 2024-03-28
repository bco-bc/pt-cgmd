/*
 * Author: André H. Juffer.
 * Created on 04/03/24, 10:09.
 *
 * Copyright (c) 2024 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/simulation.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/direction.hpp"
#include <fstream>
#include <iostream>

using namespace simploce;
using namespace simploce::param;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName =
            "/wrk3/simulation/NaClNextToChargedSurface/no-external-field-sf-vplane/NaCl-cg-surface-particle-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::obtainFrom(stream);
    stream.close();
    std::cout << "Particle specifications catalog: " << std::endl << *catalog << std::endl;

    fileName =
            "/wrk3/simulation/NaClNextToChargedSurface/no-external-field-sf-vplane/NaCl-cg-surface.ps";
    auto particleSystem = factory::particleSystem(fileName,catalog, false);
    std::cout << "Number of particles: " << particleSystem->numberOfParticles() << std::endl;

    fileName =
            "/wrk3/simulation/NaClNextToChargedSurface/no-external-field-sf-vplane/NaCl-cg-surface-interaction-parameters.dat";
    auto forceField = factory::forceField(fileName, catalog);
    std::cout << "Force field: " << std::endl << *forceField << std::endl;

    auto param = factory::simulationParameters();
    param->put<bool>("simulation.displacer.mc.in-box", true);
    param->put<bool>("simulation.forces.include-external", true);
    param->put<real_t>("simulation.forces.cutoffLR", 30.0);
    param->put<int>("simulation.nsteps", 100);
    param->put<int>("simulation.nwrite", 1);
    param->put<bool>("simulation.mesoscale", false);
    param->put<bool>("simulation.freeze-boundary", false);
    std::cout << "Simulation parameters: " << *param << std::endl;

    auto box = particleSystem->box();
    auto bc = factory::pbc_2d(box, Direction::X, Direction::Y);
    auto interactor = factory::interactor(param, forceField, bc);
    auto displacer = factory::displacer(conf::MONTE_CARLO,
                                        param,
                                        interactor,
                                        bc);

    Simulation simulation{param, particleSystem, catalog, displacer, bc, interactor};
    fileName = "/wrk3/simulation/NaClNextToChargedSurface/no-external-field-sf-vplane/simulation.dat";
    std::ofstream dataStream;
    util::open_output_file(dataStream, fileName);
    std::ofstream trajectoryStream;
    fileName = "/wrk3/simulation/NaClNextToChargedSurface/no-external-field-sf-vplane/trajectory.dat";
    util::open_output_file(trajectoryStream, fileName);
    simulation.perform(trajectoryStream, dataStream);
    dataStream.close();
    trajectoryStream.close();
}
