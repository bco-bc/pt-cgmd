/*
 * Author: André H. Juffer.
 * Created on 24/11/2021, 13:12.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "pair-potential.hpp"

#ifndef SIMULATION_NO_PAIR_POTENTIAL_HPP
#define SIMULATION_NO_PAIR_POTENTIAL_HPP

namespace simploce {

    class NoPairPotential: public pair_potential {
    public:

        NoPairPotential();

        /**
         * @return Zero energy and zero force.
         */
        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;
    };
}

#endif //SIMULATION_NO_PAIR_POTENTIAL_HPP
