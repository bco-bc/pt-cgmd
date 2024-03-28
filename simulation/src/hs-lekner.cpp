/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 7 February 2024
 */

#include "simploce/potentials/hs-lekner.hpp"
#include "simploce/potentials/lekner.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    HardSphereLekner::HardSphereLekner(ff_ptr_t forceField,
                                       bc_ptr_t bc,
                                       lekner_ptr_t lekner) :
        forceField_{std::move(forceField)}, bc_{std::move(bc)}, lekner_{std::move(lekner)}  {}

    std::pair<energy_t, force_t>
    HardSphereLekner::operator()(const simploce::p_ptr_t &pi, const simploce::p_ptr_t &pj) {
        static util::Logger logger{"simploce::HardSphereLekner::operator()"};
        logger.trace("Entering");

        static const real_t LARGE = conf::LARGE;

        // Current position.
        auto ri = pi->position();
        auto rj = pj->position();

        // Radii
        auto radius1 = pi->spec()->radius();
        auto radius2 = pj->spec()->radius();
        auto minimumDistance = radius1 + radius2;

        auto rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);
        if (  Rij <= minimumDistance() ) {
            // Hard sphere. Overlap.
            return std::move(std::make_pair<energy_t, force_t>(LARGE, force_t{LARGE, LARGE, LARGE}));
        }

        auto eps_r = forceField_->hardSphereLekner(pi->spec(), pj->spec());

        auto qi = pi->charge();
        auto qj = pj->charge();
        auto result = lekner_->forceAndEnergy(rij, qi, qj);
        result.first /= eps_r;
        result.second /= eps_r;

        logger.trace("Leaving.");
        return std::move(result);
    }


}
