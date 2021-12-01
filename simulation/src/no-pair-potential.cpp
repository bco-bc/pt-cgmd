/*
 * Author: André H. Juffer.
 * Created on 24/11/2021, 13:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/no-pair-potential.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    NoPairPotential::NoPairPotential() = default;

    std::pair<energy_t, force_t>
    NoPairPotential::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        static util::Logger logger("simploce::NoPairPotential::operator ()");
        logger.warn("No pair potential registered for ('" +
        p1->spec()->name() + "', '" + p2->spec()->name() + "') particle specification pair.");
        return std::move(std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0}));
    }
}