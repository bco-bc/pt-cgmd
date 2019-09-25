/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particle-test.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:52 PM
 */

#include "simploce/particle/ptypes.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/atomistic.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "particle-test test 1" << std::endl;
    Atomistic atomistic;
    atom_spec_ptr_t spec = ParticleSpec::createForAtom("spec_123",1.0, 2.0, 3.0);
    atom_ptr_t atom = atomistic.addAtom(1, "C", position_t{}, spec);
    std::cout << *atom << std::endl;
}

void test2() {
    std::cout << "particle-test test 2" << std::endl;
    CoarseGrained coarseGrained;
    std::cout << "particle-test test 2" << std::endl;
    bead_spec_ptr_t spec = ParticleSpec::createForBead("spec_456", 2.0, 3.0, 4.0);
    bead_ptr_t bead = coarseGrained.addBead(1, "bead1", position_t{}, spec);
    std::cout << *bead << std::endl;
    
    coarseGrained.doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
        for (auto bead : beads) {
            std::cout << *bead << std::endl;   
        }
    });
     
    //std::cout << "%TEST_FAILED% time=0 testname=test2 (particle-test) message=error message sample" << std::endl;
}

void test3() 
{
    std::cout << "particle-test test 3" << std::endl;
    CoarseGrained coarseGrained;
    bead_spec_ptr_t spec = 
            ParticleSpec::createForProtonatableBead("COOH", charge_t{-1.0}, 5.0, 3.0, 7.9);
    prot_bead_ptr_t pb = coarseGrained.addBead(1, "COOH", position_t{}, spec);
    pb->protonate();
    std::cot << "Protonated:" << std::endl;
    std::cout << "Charge: " << pb->charge() << std::endl;
    std::cout << "Charge: " << pb->charge() << std::endl;
    
    coarseGrained.doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
        for (auto bead : beads) {
            std::cout << *bead << std::endl;   
        }
    });     
}

void test4()
{
    std::cout << "particle-test test 4" << std::endl;
    CoarseGrained coarseGrained;
    bead_spec_ptr_t spec = 
            ParticleSpec::createForProtonatableBead("NH4", 0.0, 5.0, 2.0, 10.0);
    coarseGrained.addProtonatableBead(1111, "NH4", position_t{}, 1, spec);
    coarseGrained.doWithAll<void>([] (const std::vector<bead_ptr_t>& beads) {
        for (auto bead : beads) {
            std::cout << *bead << std::endl;   
        }
    });     
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% particle-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (particle-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (particle-test)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (particle-test)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (particle-test)" << std::endl;

    std::cout << "%TEST_STARTED% test3 (particle-test)\n" << std::endl;
    test3();
    std::cout << "%TEST_FINISHED% test3 (particle-test)" << std::endl;

    std::cout << "%TEST_STARTED% test4 (particle-test)\n" << std::endl;
    test4();
    std::cout << "%TEST_FINISHED% test4 (particle-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

