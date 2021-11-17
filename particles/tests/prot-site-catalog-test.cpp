/*
 * File:   prot-site-catalog-argon.cpp
 * Author: ajuffer
 *
 * Created on September 18, 2019, 3:57 PM
 */

#include "simploce/particle/protonation-site-catalog.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/p-types.hpp"
#include "simploce/util/file.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

void test1() {
    std::ifstream stream;
    
    Atomistic atomistic;

    util::open_input_file(stream, "/localdisk/resources/protonation-sites.dat");
    auto catalog = ProtonationSiteCatalog::create(stream);
    std::cout << *catalog << std::endl;
    stream.close();
    
    util::open_input_file(stream, "/localdisk/resources/particles-specs.dat");
    spec_catalog_ptr_t specs = ParticleSpecCatalog::obtainFrom(stream);
    std::cout << *specs << std::endl;
    
    spec_ptr_t specC = specs->lookup("C");
    spec_ptr_t specO = specs->lookup("O");
    spec_ptr_t specH = specs->lookup("H");
    spec_ptr_t specN = specs->lookup("N");

    
    // COOH
    auto C = atomistic.addAtom("C", specC);
    auto O1 = atomistic.addAtom("O1", specO);
    auto O2 = atomistic.addAtom("O2", specO);
    auto H = atomistic.addAtom("H", specH);
    
    // NH4
    auto N = atomistic.addAtom("N", specN);
    auto H1 = atomistic.addAtom("H1", specH);
    auto H2 = atomistic.addAtom("H2", specH);
    auto H3 = atomistic.addAtom("H3", specH);
    auto H4 = atomistic.addAtom("H4", specH);

    std::cout << "Number of atoms: " << atomistic.numberOfParticles() << std::endl;
    
    std::cout << "Looking up protonation sites...";
    /* atomistic.protonationSites(catalog);
    std::cout << "Done:" << std::endl;
    std::cout << "Number of protonation sites found: " 
              << atomistic.numberOfProtonationSites() << std::endl;    
    std::cout << "Protonation state atomistic: " << atomistic.protonationState() << std::endl;
     */
    std::cout << "Charge state atomistic: " << atomistic.charge() << std::endl;

    /*
    std::cout << "Protonating atomistic..." << std::endl;
    atomistic.doWithProtonationSites<void>([] (const std::vector<atom_prot_site_ptr_t> &sites) {
        for (const auto& site : sites) {
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
     */

    /*
    std::cout << "Deprotonating atomistic..." << std::endl;
    atomistic.doWithProtonationSites<void>([] (const std::vector<atom_prot_site_ptr_t> &sites) {
        for (const auto& site : sites) {
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
     */
}

int main(int argc, char** argv) {
    test1();
    return (EXIT_SUCCESS);
}

