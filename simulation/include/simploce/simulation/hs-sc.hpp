/*
 * Author: André H. Juffer.
 * Created on 23/11/2021, 12:50.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_HS_SC_HPP
#define SIMULATION_HS_SC_HPP

#include "pair-potential.hpp"

namespace simploce {

    /**
     * Hard sphere and screened Coulomb interaction.
     */
    class HS_SC : public pair_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param forceField Force field.
         * @param bc Boundary condition.
         * @param sc Screened Coulomb
         */
        HS_SC(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        sc_ptr_t sc_;

    };
}

#endif //SIMULATION_HS_SC_HPP
