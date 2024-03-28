/*
 * File:   leap-frog-Yiannourakou.cpp
 * Author: ajuffer
 *
 * Created on August 16, 2019, 1:00 PM
 */

#include "simploce/simulation/leap-frog.hpp"
#include "simploce/simulation/velocity-verlet.hpp"
#include "simploce/simulation/langevin-velocity-verlet.hpp"
#include "simploce/simulation/sim-data.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/types/s-types.hpp"
#include "simploce/conf/s-conf.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/param.hpp"
#include <fstream>
#include <cstdlib>
#include <iostream>

using namespace simploce;
using namespace simploce::conf;
using namespace simploce::param;

/*
 * Simple C++ Test Suite
 */

void mc(const spec_catalog_ptr_t& catalog, const ff_ptr_t& forceField) {
    std::cout << "displacer-test Monte Carlo" << std::endl;

    auto param = factory::simulationParameters();
    param->put<std::string>("simulation.displacer.mc.keep-in-box", "x,y");
    param->put<bool>("simulation.displacer.mc.in-box", false);
    param->put<bool>("simulation.mesoscale", false);
    param->put<real_t>("simulation.displacer.mc.range", 1.5);
    param->put<bool>("simulation.displacer.mc.z-non-negative", true);
    std::cout << *param << std::endl;

    auto factory = factory::particleSystemFactory(catalog);
    p_system_ptr_t particleSystem = factory->diatomic(0.12, catalog->O());
    bc_ptr_t bc = factory::pbc(particleSystem->box());
    auto interactor =
            factory::interactor(param, forceField, bc);
    auto displacer = factory::displacer(conf::MONTE_CARLO,
                                        param,
                                        interactor,
                                        bc);
    displacer->displace(particleSystem);
    displacer->displace(particleSystem);
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::obtainFrom(stream);
    stream.close();
    std::cout << *catalog << std::endl;

    fileName = "/localdisk/resources/interaction-parameters.dat";
    util::open_input_file(stream, fileName);
    auto forceField = factory::forceField(stream, catalog);
    stream.close();

    mc(catalog, forceField);

    return (EXIT_SUCCESS);
}