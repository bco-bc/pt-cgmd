/*
 * Author: André H. Juffer.
 * Created on 23/11/2021, 12:55.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/hs-sc.hpp"
#include "simploce/potentials/sc.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/simulation/s-conf.hpp"
#include <utility>

namespace simploce {

    HS_SC::HS_SC(ff_ptr_t forceField, bc_ptr_t bc) :
         forceField_{std::move(forceField)}, bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    HS_SC::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        static const real_t LARGE = conf::LARGE;

            // Current position.
        auto r1 = p1->position();
        auto r2 = p2->position();

        // Radii
        auto radius1 = p1->spec()->radius();
        auto radius2 = p2->spec()->radius();
        auto minimumDistance = radius1 + radius2;

        // Apply boundary condition.
        dist_vect_t rij = bc_->apply(r1, r2);
        auto Rij = norm<real_t>(rij);
        if ( Rij <= minimumDistance() ) {
            // Hard sphere.
            return std::move(std::make_pair<energy_t, force_t>(LARGE, force_t{LARGE, LARGE, LARGE}));
        } else {
            // Screened coulomb electrostatics.
            real_t eps_inside_rc = forceField_->hardSphereScreenedCoulomb(p1->spec(), p2->spec());
            real_t Rij2 = Rij * Rij;
            return std::move(SC::forceAndEnergy(rij,
                                                Rij,
                                                Rij2,
                                                p1->charge(),
                                                p2->charge(),
                                                eps_inside_rc));
        }

    }
}

