/*
 * Author: André H. Juffer.
 * Created on 24/11/2021, 13:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/simulation/no-pair-potential.hpp"

namespace simploce {

    NoPairPotential::NoPairPotential() = default;

    std::pair<energy_t, force_t>
    NoPairPotential::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        return std::move(std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0}));
    }
}