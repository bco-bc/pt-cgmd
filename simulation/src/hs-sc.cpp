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
#include "simploce/conf/s-conf.hpp"
#include "simploce/util/logger.hpp"
#include <utility>

namespace simploce {

    HS_SC::HS_SC(ff_ptr_t forceField, bc_ptr_t bc, std::shared_ptr<SC> screenedCoulomb) :
         forceField_{std::move(forceField)}, bc_{std::move(bc)}, screenedCoulomb_{std::move(screenedCoulomb)}  {
        util::Logger logger{"simploce::HS_SC::HS_SC()"};
        if (!forceField_) {
            util::logAndThrow(logger, "Missing force field.");
        }
        if (!bc_) {
            util::logAndThrow(logger, "Missing boundary conditions.");
        }
        if (!screenedCoulomb_) {
            util::logAndThrow(logger, "Missing screened Coulomb potential.");
        }
    }

    std::pair<energy_t, force_t>
    HS_SC::operator () (const p_ptr_t &p1, const p_ptr_t &p2) {
        util::Logger logger{"simploce::HS_SC::operator ()()"};
        logger.trace("Entering");

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
            logger.trace("Leaving");
            return std::move(std::make_pair<energy_t, force_t>(LARGE, force_t{LARGE, LARGE, LARGE}));
        } else {
            // Screened coulomb electrostatics.
            real_t eps_inside_rc = forceField_->hardSphereScreenedCoulomb(p1->spec(), p2->spec());
            real_t Rij2 = Rij * Rij;
            auto pairSC = screenedCoulomb_->forceAndEnergy(rij,
                                                           Rij,
                                                           Rij2,
                                                           p1->charge(),
                                                           p2->charge(),
                                                           eps_inside_rc);
            logger.trace("Leaving");
            return std::move(pairSC);
        }

    }
}

