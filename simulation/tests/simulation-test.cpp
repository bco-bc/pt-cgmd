/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on August 20, 2019, 12:50 PM
 */

#include "simploce/simulation/simulation.hpp"
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/interactor.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1(const spec_catalog_ptr_t& catalog, const ff_ptr_t& forceField) {
    std::cout << "simulation-test test 1" << std::endl;
    
    std::ofstream trajectoryStream, dataStream;
    util::open_output_file(trajectoryStream, "/wrk3/tests/trajectory.dat");
    util::open_output_file(dataStream, "/wrk3/tests/simulation.dat");

    auto simulationParameters = factory::simulationParameters();
    auto factory = factory::particleSystemFactory(catalog);
    auto particleSystem = factory->diatomic(0.12, catalog->O());
    auto bc = factory::boundaryCondition(particleSystem->box());
    auto interactor = factory::interactor(simulationParameters, forceField, bc);
    auto displacer = factory::displacer(conf::MONTE_CARLO, simulationParameters, interactor);
    Simulation simulation{simulationParameters, particleSystem, displacer};

    std::cout << "Simulating..." << std::endl;
    simulation.perform(trajectoryStream, dataStream);
    std::cout << "Done." << std::endl;
    
    trajectoryStream.close();
    dataStream.close();
}

void test2() {
}

int main() {
    std::string fileName = "/localdisk/resources/particles-specs.dat";
    auto catalog = factory::particleSpecCatalog(fileName);
    fileName = "/localdisk/resources/interaction-parameters.dat";
    auto forceField = factory::forceField(fileName, catalog);
    test1(catalog, forceField);

    return (EXIT_SUCCESS);
}

