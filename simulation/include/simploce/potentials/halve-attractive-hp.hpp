/**
 * Author: André H. Juffer.
 * Created on 14/11/2023
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */


#ifndef SIMULATION_HALVE_ATTRACTIVE_HP_HPP
#define SIMULATION_HALVE_ATTRACTIVE_HP_HPP

#include "pair-potential.hpp"
#include "../types/s-types.hpp"

namespace simploce {

    /**
     * Halve attractive quartic bonded potential, U(r) = 0.5 * k * (r - r0)^2, where r is a distance,
     * k is a force constant (fc) and r0 is the equilibrium distance. Halve attractive implies
     * U(r) = 0 and dU(r)/dr = 0 for r <= r0.
     */
    class HalveAttractiveHP : public pair_potential {
    public:
        HalveAttractiveHP(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
    };
}

#endif //SIMULATION_HALVE_ATTRACTIVE_HP_HPP
