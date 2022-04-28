/*
 * Author: André H. Juffer.
 * Created on 01/12/2021, 12:28.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_LJ_SF_HPP
#define SIMULATION_LJ_SF_HPP

#include "pair-potential.hpp"
#include "simploce/types/s-types.hpp"

namespace simploce {

    /**
     * Lennard-Jones plus shifted force electrostatics.
     */
    class LJ_SF: public pair_potential {
    public:

        LJ_SF(ff_ptr_t forceField, bc_ptr_t bc, sf_ptr_t sf);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        sf_ptr_t sf_;
    };
}

#endif //SIMULATION_LJ_SF_HPP
