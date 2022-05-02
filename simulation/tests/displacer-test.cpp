/*
 * File:   leap-frog-test.cpp
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
using namespace simploce::param;

/*
 * Simple C++ Test Suite
 */

void test1(const spec_catalog_ptr_t& catalog, const ff_ptr_t& forceField) {
    std::cout << "displacer-test test 1" << std::endl;

    auto simulationParameters = factory::simulationParameters();
    std::cout << *simulationParameters << std::endl;
    
    box_ptr_t box = factory::box(length_t{5.0});
    bc_ptr_t bc = factory::boundaryCondition(box);

    auto factory = factory::particleSystemFactory(catalog);
    p_system_ptr_t particleSystem = factory->diatomic(0.12, catalog->O());
    auto interactor =
            factory::interactor(simulationParameters, forceField, bc);
    auto displacer = factory::displacer(conf::MONTE_CARLO, simulationParameters, interactor, bc);
    auto result = interactor->interact(particleSystem);
    std::cout << "BEFORE: Non-bonded: " << std::get<1>(result);
    std::cout << ", bonded: " << std::get<0>(result) << std::endl;
    std::cout << "; external: " << std::get<2>(result) << std::endl;
    std::cout << "Displacing particle system..." << std::endl;
    for (int k = 0; k != 100; ++k) {
        displacer->displace(particleSystem);
    }
    std::cout << "Done." << std::endl;
    result = interactor->interact(particleSystem);
    std::cout << "AFTER: Non-bonded: " << std::get<1>(result);
    std::cout << ", bonded: " << std::get<0>(result) << std::endl;
    std::cout << "; external: " << std::get<2>(result) << std::endl;
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

    test1(catalog, forceField);

    return (EXIT_SUCCESS);
}