/*
 * File:   protonation-site-argon.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 18, 2019, 3:55 PM
 */

#include "simploce/particle/protonation-site.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/p-types.hpp"
#include "simploce/particle/atomistic.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    Atomistic atomistic;
    
    // Deprotonated specification
    spec_ptr_t O_dp =
            ParticleSpec::create("O", -1.0, 1.0, 1.0, false, "O deprotonated");
    spec_ptr_t O_p =
            ParticleSpec::create("O", -0.5, 1.0, 1.0, false, "O protonated");

    // Protonated specifications.
    spec_ptr_t H_dp =
            ParticleSpec::create("H", 0.0, 1.0, 1.0, false, "H deprotonated");
    spec_ptr_t H_p =
            ParticleSpec::create("H", 0.5, 1.0, 1.0, false, "H protonated");

    // Atoms.
    atom_ptr_t atom1 = atomistic.addAtom("O", O_p);
    atom_ptr_t atom2 = atomistic.addAtom("H", H_p);

    std::vector<atom_ptr_t> atoms{atom1, atom2};
    std::vector<spec_ptr_t> deprotonated{O_dp, H_dp};
    std::vector<spec_ptr_t> protonated{O_p, H_p};

    std::vector<id_pair_t> bonds{};
    id_pair_t pair = std::make_pair(atom1->id(), atom2->id());
    bonds.push_back(pair);
    ProtonationSite<Atom> site("argon", atoms, bonds, deprotonated, protonated);
    
    std::cout << "Is protonated: " << site.isProtonated() << std::endl;
    std::cout << "Total charge: " << site.charge() << std::endl;
    std::cout << "Protonation state: " << site.protonationState() << std::endl;
    std::cout << std::endl;

    if ( !site.isProtonated() ) {
        site.protonate();
    }
    std::cout << "Is protonated: " << site.isProtonated() << std::endl;
    std::cout << "Total charge: " << site.charge() << std::endl;
    std::cout << "Protonation state: " << site.protonationState() << std::endl;
    std::cout << std::endl;

    if (site.isProtonated() ) {
        site.deprotonate();
    }

    std::cout << "Is protonated: " << site.isProtonated() << std::endl;
    std::cout << "Total charge: " << site.charge() << std::endl;
    std::cout << "Protonation state: " << site.protonationState() << std::endl;
}

int main() {
    test1();
    return (EXIT_SUCCESS);
}

