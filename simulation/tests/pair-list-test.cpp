/*
 * File:   pair-list-test.cpp
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 4, 2019, 2:55 PM
 */

#include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/simulation/cell-pair-list-generator.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/param.hpp"
#include <cstdlib>
#include <iostream>
#include <string>
#include <chrono>

using namespace simploce;
using namespace simploce::param;
using namespace std::chrono;

void test2(const param_ptr_t& param, p_system_ptr_t &coarseGrained) {
    auto box = coarseGrained->box();
    auto bc = factory::boundaryCondition(box);
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(param, bc);
    auto pairLists = generator->generate(coarseGrained);
    const auto& particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

void distanceBased(const param_ptr_t& param, p_system_ptr_t & particleSystem) {
    auto box = particleSystem->box();
    auto bc = factory::boundaryCondition(box);
    DistancePairListGenerator generator(param, bc);
    auto pairLists = generator.generate(particleSystem);
    const auto& particlePairs = pairLists.particlePairList();
    std::cout << "Distance-based number of pairs: " << particlePairs.size() << std::endl;
}

void cellBased(const param_ptr_t& param, p_system_ptr_t& particleSystem) {
    auto box = particleSystem->box();
    auto bc = factory::boundaryCondition(box);
    CellPairListGenerator generator(param, bc);
    auto pairLists = generator.generate(particleSystem);
    const auto& particlePairs = pairLists.particlePairList();
    std::cout << "Cell-based number of pairs: " << particlePairs.size() << std::endl;
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGINFO);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    auto catalog = ParticleSpecCatalog::obtainFrom(stream);
    stream.close();
    std::clog << "Particle specs: " << *catalog << std::endl;

    auto factory = factory::protonatableParticleSystemFactory(catalog);
    auto param = factory::simulationParameters();

    //auto polarizableWater = factory->polarizableWater(factory::box(7.27));
    // test2(param, polarizableWater);

    //auto argon = factory->argon(factory::box(3.47786));
    //distanceBased(argon);

    //auto electrolyte = factory->simpleElectrolyte();
    //auto electrolyte= factory::particleSystem("/home/andre/simulations/electrolyte/electrolyte-mc-3.ps",
    //                                          catalog,
    //                                          false);
    //distanceBased(param, electrolyte);



    // WARNING: Mesoscale units!


    /*
    auto box = factory::box(10.0, 10., 40.0);
    auto particleSystem =
            factory->polymerSolution(box,
                                     15,
                                     "PMU",
                                     400,
                                     1.0,
                                     6000,
                                     "H2Om",
                                     2.0,
                                     true);
    */
    /*
    auto box = factory::box(10, 10, 20);
    auto particleSystem = factory->dropletPolymerSolution(box,
                                            15,
                                            "PMU",
                                            0,
                                            1.0,
                                            3000,
                                            "H2Om",
                                            3000,
                                            "DROm",
                                            1.0);
                                            */

    std::cout << "Reading particle system..." << std::endl;
    auto particleSystem =
            simploce::factory::particleSystem("/wrk3/tests/droplets-polymer-solution.ps",
                                              catalog,
                                              true);
    std::cout << "Number of particles: " << particleSystem->numberOfParticles() << std::endl;

    param->put("simulation.timestep", 0.01);
    param->put("simulation.mesoscale", true);
    param->put("simulation.forces.cutoff", 1.0);
    param->put("simulation.temperature", 1.0);

    param::write(std::clog, *param);

    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::cout << std::endl;
    std::cout << "Cell-based (with setup):" << std::endl;
    auto start = high_resolution_clock::now();
    cellBased(param, particleSystem);
    auto stop = high_resolution_clock::now();
    auto durationCellBased = duration_cast<microseconds>(stop - start);
    std::cout << "Duration: " << durationCellBased.count() << std::endl;

    std::cout << std::endl;
    std::cout << "Distance-based:" << std::endl;
    start = high_resolution_clock::now();
    distanceBased(param, particleSystem);
    stop = high_resolution_clock::now();
    auto durationDistanceBased = duration_cast<microseconds>(stop - start);
    std::cout << "Duration: " << durationDistanceBased.count() << std::endl;

    std::cout << std::endl;
    std::cout << "Cell-based (no setup):" << std::endl;
    start = high_resolution_clock::now();
    cellBased(param, particleSystem);
    stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Duration: " << duration.count() << std::endl;
    std::cout << std::endl;

    std::cout << "Ratio durations cell-based (with setup) / distance_based: ";
    std::cout  << real_t(durationCellBased.count()) / real_t(durationDistanceBased.count()) << std::endl;

    std::cout << "Ratio durations cell-based (no setup) / distance_based: ";
    std::cout  << real_t(duration.count()) / real_t(durationDistanceBased.count()) << std::endl;

    return (EXIT_SUCCESS);
}

