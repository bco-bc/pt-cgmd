/*
 * File:   protonation-site-test.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:55 PM
 */

#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/ptypes.hpp"
#include "simploce/particle/atomistic.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "protonation-site-test test 1" << std::endl;
    
    Atomistic atomistic;
    
    // Identical spec names!
    spec_ptr_t spec1_dp = ParticleSpec::create("pspec1", 0.0, 1.0, 1.0);
    spec_ptr_t spec1_p = ParticleSpec::create("pspec1", 1.0, 1.0, 1.0);
    
    atom_ptr_t atom1 = atomistic.addAtom(1, "test1", position_t{}, spec1_dp);
    atom_ptr_t atom2 = atomistic.addAtom(11, "test2", position_t{}, spec1_dp);
    
    std::vector<atom_ptr_t> atoms{atom1, atom2};    
    std::vector<spec_ptr_t> deprotonated{spec1_dp, spec1_dp};    
    std::vector<spec_ptr_t> protonated{spec1_p, spec1_p};
    
    ProtonationSite<Atom> site("test", atoms, deprotonated, protonated);
    
    std::cout << "Is protonated: " << site.isProtonated() << std::endl;
    
    if ( !site.isProtonated() ) {
        site.protonate();
    }
    std::cout << "Is protonated: " << site.isProtonated() << std::endl;
    if ( site.isProtonated() ) {
        site.deprotonate();
    }
    std::cout << "Is protonated: " << site.isProtonated() << std::endl;
    std::cout << "Protonation state: " << site.protonationState() << std::endl;
}

void test2() {
    std::cout << "protonation-site-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (protonation-site-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% protonation-site-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (protonation-site-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (protonation-site-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test2 (protonation-site-test)\n" << std::endl;
    //test2();
    //std::cout << "%TEST_FINISHED% time=0 test2 (protonation-site-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

