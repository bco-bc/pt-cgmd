/*
 * File:   ext-potential-test.cpp
 * Author: André H. Juffer, Biocenter Oulu
 *
 * Created on October 16, 2019, 4:08 PM
 */

#include "simploce/forcefields/ext-potential.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/atomistic.hpp"
#include <iostream>

using namespace simploce;

class TestExtPotential : public ext_potential<Atom> {
public:

    using particle_ptr_t = ext_potential<Atom>::particle_ptr_t;

    energy_t energy(const particle_ptr_t& p) const override {
        return 1.0;
    }

    energy_t force(particle_ptr_t& p) const override {
        return 0.0;
    }
};

int main() {
    Atomistic atomistic;

    auto spec = ParticleSpec::create("test", 0.0, 1.0, 0.1, "test");
    atomistic.addAtom(1, "test", position_t{}, spec);
    atomistic.addAtom(2, "test", position_t{}, spec);

    TestExtPotential potential;
    auto energy = atomistic.doWithAll<energy_t>([potential] (const std::vector<atom_ptr_t> & atoms) {
        energy_t energy = 0.0;
        for (const auto& atom : atoms) {
            energy += potential.energy(atom);
        }
        return energy;
    });
    std::cout << "Number of atoms: " << atomistic.numberOfParticles() << std::endl;
    std::cout << "Total energy: " << energy << std::endl;
    return 0;
}