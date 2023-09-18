/*
 * Author: André H. Juffer.
 * Created on 24/11/2021, 13:12.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/pair-potential.hpp"

#ifndef SIMULATION_NON_INTERACTING_HPP
#define SIMULATION_NON_INTERACTING_HPP

namespace simploce {

    /**
     * Non-interacting pair potential. That is, interaction energy and forces vector components are zero.
     */
    class NonInteracting: public pair_potential {
    public:

        NonInteracting();

        /**
         * @return Zero energy and zero force.
         */
        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;
    };
}

#endif //SIMULATION_NON_INTERACTING_HPP
