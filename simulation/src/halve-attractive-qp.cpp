/*
 * Author: André H. Juffer.
 * Created on 16/11/2021, 18:49.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/halve-attractive-qp.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"


namespace simploce {

    HalveAttractiveQP::HalveAttractiveQP(ff_ptr_t forceField, bc_ptr_t bc) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    HalveAttractiveQP::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        // Get r0 and fc.
        auto params = forceField_->halveAttractiveQuartic(p1->spec(), p2->spec());
        auto r0 = params.first;
        auto fc = params.second;

        // Current positions.
        const auto &r1 = p1->position();
        const auto &r2 = p2->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        real_t dR = Rij - r0;

        if ( dR > 0.0 ) {
            real_t dR3 = dR * dR * dR;
            real_t dR4 = dR * dR3;
            energy_t energy{0.5 * fc * dR4};
            real_t derQP = 2.0 * fc * dR3;
            dist_vect_t unitVector = rij / Rij;
            force_t f{};
            for (std::size_t k = 0; k != 3; ++k) {
                f[k] = -derQP * unitVector[k];
            }
            return  std::move(std::make_pair(energy, f));
        } else {
            return std::move(std::make_pair(energy_t{0.0}, force_t{0.0, 0.0, 0.0}));
        }
    }

}

