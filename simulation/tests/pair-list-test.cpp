/*
 * File:   pair-list-test.cpp
 * Author: ajuffer
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
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

void test2(p_system_ptr_t &coarseGrained) {
    auto box = coarseGrained->box();
    auto bc = factory::boundaryCondition(box);
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(bc);
    auto pairLists = generator->generate(coarseGrained);
    auto particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

void test3(p_system_ptr_t & atomistic) {
    auto box = atomistic->box();
    auto bc = factory::boundaryCondition(box);
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(bc);
    auto pairLists = generator->generate(atomistic);
    auto particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    auto catalog = ParticleSpecCatalog::obtainFrom(stream);

    auto factory = factory::protonatableParticleSystemFactory(catalog);

    auto polarizableWater = factory->polarizableWater(factory::box(7.27));
    test2(polarizableWater);

    //auto argon = factory->argon(factory::box(3.47786));
    //test3(argon);

    return (EXIT_SUCCESS);
}

