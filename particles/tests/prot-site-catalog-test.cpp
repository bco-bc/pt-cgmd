/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   prot-site-catalog-test.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:57 PM
 */

#include "simploce/particle/protonation-site-catalog.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/ptypes.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "prot-site-catalog-test test 1" << std::endl;
    std::ifstream stream;
    
    Atomistic atomistic;

    file::open_input(stream, "/home/ajuffer/simploce/particles/resources/protonation-sites.dat");
    prot_site_catalog_ptr_t catalog = ProtonationSiteCatalog::create(stream);
    std::cout << *catalog << std::endl;
    stream.close();
    
    file::open_input(stream, "/home/ajuffer/simploce/particles/resources/particles-specs.dat");
    spec_catalog_ptr_t specs = ParticleSpecCatalog::create(stream);
    std::cout << *specs << std::endl;
    
    atom_spec_ptr_t specC = specs->lookup("C");
    atom_spec_ptr_t specO = specs->lookup("O");
    atom_spec_ptr_t specH = specs->lookup("H");
    atom_spec_ptr_t specN = specs->lookup("N");
   
    position_t r;
    
    // COOH
    atom_ptr_t C = atomistic.addAtom("C", r, specC);
    atom_ptr_t O1 = atomistic.addAtom("O1", r, specO);
    atom_ptr_t O2 = atomistic.addAtom("O2", r, specO);
    atom_ptr_t H = atomistic.addAtom("H", r, specH);
    
    // NH4
    atom_ptr_t N = atomistic.addAtom("N", r, specN);
    atom_ptr_t H1 = atomistic.addAtom("H1", r, specH);
    atom_ptr_t H2 = atomistic.addAtom("H2", r, specH);
    atom_ptr_t H3 = atomistic.addAtom("H3", r, specH);
    atom_ptr_t H4 = atomistic.addAtom("H4", r, specH);

    std::cout << "Number of atoms: " << atomistic.size() << std::endl;
    
    std::cout << "Looking up protonation sites...";
    atomistic.protonationSites(catalog);
    std::cout << "Done:" << std::endl;
    std::cout << "Number of protonation sites found: " 
              << atomistic.numberOfProtonationSites() << std::endl;    
    std::cout << "Protonation state atomistic: " << atomistic.protonationState() << std::endl;
    std::cout << "Charge state atomistic: " << atomistic.charge() << std::endl;
    
    std::cout << "Protonating atomistic..." << std::endl;
    atomistic.operateOnProtonationSites<void>([] (const std::vector<atom_prot_site_ptr_t> &sites) {
        for (auto site : sites) {
            std::cout << "Site (before): " << site->name() 
                      << ", protonation state: " << site->protonationState()
                      << ", charge state: " << site->charge() << std::endl;
            if ( !site->isProtonated() ) {
                site->protonate();
            }
            std::cout << "Site (after): " << site->name() 
                      << ", protonation state: " << site->protonationState()
                      << ", charge state: " << site->charge() << std::endl;
        }
    });
    std::cout << "Done:" << std::endl;
    std::cout << "Protonation state atomistic: " << atomistic.protonationState() << std::endl;
    std::cout << "Charge state atomistic: " << atomistic.charge() << std::endl;

    std::cout << "Deprotonating atomistic..." << std::endl;
    atomistic.operateOnProtonationSites<void>([] (const std::vector<atom_prot_site_ptr_t> &sites) {
        for (auto site : sites) {
            std::cout << "Site (before): " << site->name()
                      << ", protonation state: " << site->protonationState() 
                      << ", charge state: " << site->charge() << std::endl;
            if ( site->isProtonated() ) {
                site->deprotonate();
            }
            std::cout << "Site (after): " << site->name() 
                      << ", protonation state: " << site->protonationState()
                      << ", charge state: " << site->charge() << std::endl;
        }
    });
    std::cout << "Done:" << std::endl;
    std::cout << "Protonation state atomistic: " << atomistic.protonationState() << std::endl;
    std::cout << "Charge state atomistic: " << atomistic.charge() << std::endl;
}

void test2() {
    std::cout << "prot-site-catalog-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (prot-site-catalog-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% prot-site-catalog-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (prot-site-catalog-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (prot-site-catalog-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test2 (prot-site-catalog-test)\n" << std::endl;
    //test2();
    //std::cout << "%TEST_FINISHED% time=0 test2 (prot-site-catalog-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

