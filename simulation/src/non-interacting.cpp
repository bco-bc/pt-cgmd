/*
 * Author: André H. Juffer.
 * Created on 24/11/2021, 13:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/non-interacting.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    NonInteracting::NonInteracting() = default;

    std::pair<energy_t, force_t>
    NonInteracting::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        static energy_t energy{0.0};
        static force_t f{0.0, 0.0, 0.0};
        return std::move(std::make_pair(energy, f));
    }
}