/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 26, 2021., 11:42 PM.
 */

#include "simploce/forcefields/static-uniform.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bead.hpp"
#include <cstdlib>

using namespace simploce;

void test1() {
    CoarseGrained cg;

    auto spec = ParticleSpec::create("test", 0.2, 1.0, 0.1, "test");
    auto bead = cg.addBead(1, "test", position_t{}, spec);

    srf_charge_density_t sigma{0.1};
    rel_perm_t eps_r{4.0};
    StaticUniform<Bead> potential{sigma, eps_r};
    std::cout << "Surface charge density: " << sigma << std::endl;
    std::cout << "Relative permittivity: " << eps_r << std::endl;
    std::cout << "Bead charge value: " << bead->charge() << std::endl;
    bead->position(position_t{1.0, 0.0, 0.0});
    std::cout << "Bead: " << *bead << std::endl;
    std::cout << "Potential at x=1.0: " << potential.energy(bead) << std::endl;
    potential.force(bead);
    std::cout << "Force at x=1.0: " <<bead->force() << std::endl;
    bead->position(position_t{-0.5, 0.0, 0.0});
    std::cout << "Bead: " << *bead << std::endl;
    std::cout << "Potential at x=-0.5: " << potential.energy(bead) << std::endl;
    potential.force(bead);
    std::cout << "Force at x=-0.5: " << bead->force() << std::endl;
}

int main() {
    test1();
    return (EXIT_SUCCESS);
}

