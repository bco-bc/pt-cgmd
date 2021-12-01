/*
 * Author: André H. Juffer.
 * Created on 13/11/2021, 17:23.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_HP_HPP
#define SIMULATION_HP_HPP

#include "simploce/potentials/pair-potential.hpp"

namespace simploce {

    /**
     * Harmonic potential, U(r) = 0.5 * k * (r - r0)^2, where r is a distance, k is a
     * force constant (fc) and r0 is the equilibrium distance.
     */
    class HP : public pair_potential {
    public:

        HP(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &pi, const p_ptr_t &pj) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
    };

}

#endif //SIMULATION_HP_HPP
