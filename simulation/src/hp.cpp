/*
 * Author: André H. Juffer.
 * Created on 13/11/2021, 17:23.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/hp.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"

namespace simploce {

    HP::HP(ff_ptr_t forceField, bc_ptr_t bc) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    HP::operator () (const p_ptr_t &pi, const p_ptr_t &pj) {
        // Get r0 and fc parameters.
        auto params = forceField_->harmonic(pi->spec(), pj->spec());
        auto r0 = params.first;
        auto fc = params.second;

        // Current positions.
        const auto &ri = pi->position();
        const auto &rj = pj->position();

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);
        real_t dR = Rij - r0;

        // Potential energy, kJ/mol
        energy_t energy{0.5 * fc * dR * dR};

        // Force on particle #1.
        dist_vect_t unitVector = rij / Rij;
        real_t dHPdR = fc * dR;
        force_t f{};
        for (std::size_t k = 0; k != 3; ++k) {
            f[k] = -dHPdR * unitVector[k];
        }

        // Done.
        return std::move(std::make_pair(energy, f));
    }

}
