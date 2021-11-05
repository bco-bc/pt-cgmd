/*
 * Author: André H. Juffer.
 * Created on 27/10/2021, 15:06.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/forcefields/lj.hpp"
#include "simploce/forcefields/lj-electrostatic.hpp"
#include "simploce/forcefields/no-bc.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-spec.hpp"
#include <cstdlib>

using namespace simploce;

void test1(CoarseGrained& cg) {
    using param_t = LennardJones<Bead>::param_t;
    param_t param;
    param.add("test", "test", std::make_pair(1.298e-03, 1.298e-03));
    LennardJones<Bead,no_bc> lj(param);
    auto energy = cg.doWithAll<energy_t>([lj] (const std::vector<bead_ptr_t>& beads) {
        auto bead1 = beads[0];
        auto bead2 = beads[1];
        energy_t energy = lj.energy(bead1, bead2);
        //energy = lj.force(bead1, bead1);
        return std::move(energy);
    });
    std::cout << "LJ interaction energy: " << energy << std::endl;
}

void test2(CoarseGrained& cg) {
    LennardJonesElectrostatic<Bead> ljElec;
    auto energy = cg.doWithAll<energy_t>([ljElec] (const std::vector<bead_ptr_t>& beads) {
        auto bead1 = beads[0];
        auto bead2 = beads[1];
        energy_t energy = ljElec.energy(bead1, bead2);
        energy = ljElec.force(bead1, bead1);
        return std::move(energy);
    });
    std::cout << "LJ+electrostatic interaction energy: " << energy << std::endl;
}

int main() {
    CoarseGrained cg;
    auto spec1 = ParticleSpec::create("test", 1.0, 1.0, 0.1, "test");
    auto spec2 = ParticleSpec::create("test", 1.0, 1.0, 0.1, "test");

    cg.addBead(1, "p1", position_t{0.0, 0.0, 0.0}, spec1);
    cg.addBead(2, "p2", position_t{0.3, 0.0, 0.0}, spec2);
    test1(cg);
    test2(cg);
    return (EXIT_SUCCESS);
}