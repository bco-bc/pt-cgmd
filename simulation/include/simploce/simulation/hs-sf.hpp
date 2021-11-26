/*
 * Author: André H. Juffer.
 * Created on 22/11/2021, 21:39.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_HS_SF_HPP
#define SIMULATION_HS_SF_HPP

#include "pair-potential.hpp"

namespace simploce {

    /**
     * Hard sphere potential, Shifted force Coulomb interaction. Should not be used when
     * force calculations are required.
     */
    class HS_SF: public pair_potential {
    public:

        /**
         * Constructor. All argument are required.
         * @param forceField Force field.
         * @param bc Boundary condition.
         * @param sf Shifted force electrostatic interaction.
         */
        HS_SF(ff_ptr_t forceField, bc_ptr_t bc, sf_ptr_t sf);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        sf_ptr_t sf_;
    };
}

#endif //SIMULATION_HS_SF_HPP
