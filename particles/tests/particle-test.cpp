/*
 * File:   particle-argon.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 18, 2019, 3:52 PM
 */

#include "simploce/particle/p-types.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/atomistic.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void testContinuous() {
    std::cout << "particle-test argon 1" << std::endl;
    Atomistic atomistic;
    auto spec = ParticleSpec::create("spec_123",1.0, 2.0, 3.0, "argon");
    auto atom = atomistic.addAtom("C", spec);
    force_t f{1.0, 2.0, 3.0};
    atom->force(f);
    std::cout << *atom << std::endl;
    std::cout << "Force: " << atom->force() << std::endl;
    atom->resetForce();
    std::cout << "Force after reset: " << atom->force() << std::endl;
}

void test2() {
    std::cout << "particle-test argon 2" << std::endl;
    CoarseGrained coarseGrained;
    std::cout << "particle-test argon 2" << std::endl;
    auto spec = ParticleSpec::create("spec_456", 2.0, 3.0, 4.0, "test2");
    auto bead = coarseGrained.addBead("bead1", spec);
    std::cout << *bead << std::endl;
    
    coarseGrained.doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
        for (const auto& bead : beads) {
            std::cout << *bead << std::endl;   
        }
    });
}

void test3() 
{
    std::cout << "particle-test argon 3" << std::endl;
    CoarseGrained coarseGrained;
    auto spec = ParticleSpec::create("COOH", charge_t{-1.0}, 5.0, 3.0, 7.9, "test3");
    auto bead =
        coarseGrained.addBead("COOH", spec);
    std::cout << "Charge: " << bead->charge() << std::endl;
    std::cout << "Charge: " << bead->charge() << std::endl;
    
    coarseGrained.doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
        for (const auto& bead : beads) {
            std::cout << *bead << std::endl;   
        }
    });     
}

void test4()
{
    std::cout << "particle-argon test 4" << std::endl;
    CoarseGrained coarseGrained;
    auto spec =
        ParticleSpec::create("NH4", 0.0, 5.0, 2.0, 10.0, "argon");
    coarseGrained.addBead("NH4", spec);
    coarseGrained.doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
        for (const auto& bead : beads) {
            std::cout << *bead << std::endl;
            std::cout << "Is ion? " << bead->isIon() << std::endl;
        }
    });
}

int main() {
    testContinuous();
    test2();
    test3();
    test4();
    return (EXIT_SUCCESS);
}

