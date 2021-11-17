/*
 * File:   pair-list-test.cpp
 * Author: ajuffer
 *
 * Created on September 4, 2019, 2:55 PM
 */

// #include "simploce/simulation/cell-lists.hpp"
// #include "simploce/simulation/distance-lists.hpp"
#include "simploce/simulation/distance-pair-list-generator.hpp"
#include "simploce/simulation/pair-lists.hpp"
#include "simploce/simulation/s-factory.hpp"
// #include "simploce/simulation/s-types.hpp"
// #include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/simulation/protonatable-particle-system-factory.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/simulation/s-properties.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace simploce;

/*
void test1(const spec_catalog_ptr_t &catalog) {
    std::cout << "pair-list-test test 1" << std::endl;
    
    using cell_generator_t = CellLists<Bead>;   
    using dist_generator_t = DistanceLists<Bead>; 
    using pairlists_t = PairLists<Bead>;
        

    auto box = factory::box(7.27);
    auto smf = factory::simulationModelFactory(catalog);
    auto sm = smf->polarizableWater(box, 997.0479, 298.15, 2560);
    
    std::clog << "Number of beads: " << sm->size() << std::endl;
    std::clog << std::endl;
    
    dist_generator_t distGen(sm->box(), sm->boundaryCondition());

    auto start = std::chrono::system_clock::now();
    auto pairlists = 
        sm->doWithAllFreeGroups<pairlists_t>([distGen] (const std::vector<bead_ptr_t>& all,
                                                        const std::vector<bead_ptr_t>& free,
                                                        const std::vector<bead_group_ptr_t>& groups) {
        return distGen.generate(all, free, groups);
    });
    auto end = std::chrono::system_clock::now();
    auto distDiff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::clog << "Distance-based pair lists: Time used: " << distDiff << std::endl;   
    std::clog << std::endl;

    cell_generator_t cellGen(sm->box(), sm->boundaryCondition());
    //cellGen.sideLength(0.25 * util::cutoffDistance(box));
    start = std::chrono::system_clock::now();
    pairlists = 
        sm->doWithAllFreeGroups<pairlists_t>([cellGen] (const std::vector<bead_ptr_t>& all,
                                                        const std::vector<bead_ptr_t>& free,
                                                        const std::vector<bead_group_ptr_t>& groups) {
        return cellGen.generate(all, free, groups);
    });
    end = std::chrono::system_clock::now();
    auto cellDiff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::clog << "Cell-based pair lists: Time used: " << cellDiff << std::endl;  
    std::clog << "Speed ratio (distance/cell): " << real_t(distDiff)/real_t(cellDiff)
              << std::endl;
    std::clog << std::endl;
}
*/

void test2(cg_sys_ptr_t &coarseGrained) {
    auto box = coarseGrained->box();
    auto bc = factory::pbc(box);
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(box, bc);
    auto pairLists =
            coarseGrained->doWithAllFreeGroups<PairLists>([generator] (
                    const std::vector<p_ptr_t>& all,
                    const std::vector<p_ptr_t>& free,
                    const std::vector<pg_ptr_t>& groups) {
        return generator->generate(all, free, groups);
    });
    auto particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

void test3(at_sys_ptr_t & atomistic) {
    auto box = atomistic->box();
    auto bc = factory::pbc(box);
    pair_list_gen_ptr_t generator = factory::pairListsGenerator(box, bc);
    auto pairLists =
            atomistic->doWithAllFreeGroups<PairLists>([generator] (
                    const std::vector<p_ptr_t>& all,
                    const std::vector<p_ptr_t>& free,
                    const std::vector<pg_ptr_t>& groups) {
        return generator->generate(all, free, groups);
    });
    auto particlePairs = pairLists.particlePairList();
    std::cout << "Number of pairs: " << particlePairs.size() << std::endl;
}

int main(int argc, char** argv) {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    auto catalog = ParticleSpecCatalog::obtainFrom(stream);
    //std::clog << *catalog << std::endl;

    //test1(catalog);

    auto factory = factory::protonatableParticleSystemFactory(catalog);
    auto polarizableWater = factory->polarizableWater(factory::box(7.27));
    test2(polarizableWater);

    auto argon = factory->argon(factory::box(3.47786));
    test3(argon);

    return (EXIT_SUCCESS);
}

