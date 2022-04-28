/*
 * File:   pair-list-test.cpp
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 4, 2019, 2:55 PM
 */

#include "simploce/simulation/distance-pair-list-generator.hpp"
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

using namespace simploce;
using namespace simploce::param;

void test2(const param_ptr_t& param, p_system_ptr_t &coarseGrained) {
    auto box = coarseGrained->box();
    auto bc = factory::boundaryCondition(box);
    auto cutoff = param->get<real_t>("forces.nb.cutoff");
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(cutoff, bc);
    auto pairLists = generator->generate(coarseGrained);
    const auto& particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

void test3(const param_ptr_t& param, p_system_ptr_t & atomistic) {
    auto box = atomistic->box();
    auto bc = factory::boundaryCondition(box);
    auto cutoff = param->get<real_t>("forces.nb.cutoff");
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(cutoff, bc);
    auto pairLists = generator->generate(atomistic);
    const auto& particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    auto catalog = ParticleSpecCatalog::obtainFrom(stream);

    auto factory = factory::protonatableParticleSystemFactory(catalog);
    auto param = factory::simulationParameters();
    param->put<real_t>("forces.nb.cutoff", 2.0);

    auto polarizableWater = factory->polarizableWater(factory::box(7.27));
    test2(param, polarizableWater);
    std::cout << "Parameters: " << std::endl << *param << std::endl;

    //auto argon = factory->argon(factory::box(3.47786));
    //test3(argon);

    //auto electrolyte = factory->simpleElectrolyte();
    //auto electrolyte= factory::particleSystem("/home/andre/simulations/electrolyte/electrolyte-mc-3.ps",
    //                                          catalog,
    //                                          false);
    //test3(param, electrolyte);

    return (EXIT_SUCCESS);
}

