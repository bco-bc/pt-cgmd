/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 6/1/2022.
 */

#include "simploce/potentials/hp-sr.hpp"
#include "simploce/potentials/soft-repulsion.hpp"
#include "simploce/potentials/hp.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/logger.hpp"
#include <stdexcept>
#include <limits>

namespace simploce {

    HarmonicSoftRepulsion::HarmonicSoftRepulsion(ff_ptr_t forceField,
                                                 bc_ptr_t bc,
                                                 std::shared_ptr<SoftRepulsion> softRepulsion,
                                                 std::shared_ptr<HP> harmonic)
        : forceField_{std::move(forceField)}, bc_{std::move(bc)},
          softRepulsion_{std::move(softRepulsion)}, harmonic_{std::move(harmonic)} {
        util::Logger logger{"simploce::HarmonicSoftRepulsion::HarmonicSoftRepulsion()"};
        if (!forceField_) {
            util::logAndThrow(logger, "Missing force field.");
        }
        if (!bc_) {
            util::logAndThrow(logger, "Missing boundary condition.");
        }
        if (!softRepulsion_) {
            util::logAndThrow(logger, "Missing soft repulsion potential.");
        }
        if (!harmonic) {
            util::logAndThrow(logger, "Missing harmonic potential.");
        }
    }

    std::pair<energy_t, force_t>
    HarmonicSoftRepulsion::operator()(const p_ptr_t &p1, const p_ptr_t &p2) {
        static util::Logger logger("simploce::HarmonicSoftRepulsion::operator()()");
        logger.trace("Entering.");

        logger.debug(p1->name() + ", " + p2->name() + ": Particle pair.");

        static auto eps = real_t(std::numeric_limits<float>::min());
        logger.debug(std::to_string(eps) + ": Minimum value for distance between two particles.");

        // Get interaction parameters.
        auto params = forceField_->harmonicSoftRepulsion(p1->spec(), p2->spec());
        auto r0 = std::get<0>(params);
        auto fc = std::get<1>(params);
        auto Aij = std::get<2>(params);

        // Current positions.
        const auto &ri = p1->position();
        const auto &rj = p2->position();

        // Distance.
        dist_vect_t rij = bc_->apply(ri, rj);
        auto Rij = norm<real_t>(rij);
        auto Rij2 = Rij * Rij;

        // Soft repulsion.
        auto pairSR = softRepulsion_->forceAndEnergy(rij, Rij, Rij2, Aij);

        // Harmonic.
        auto pairHP = harmonic_->forceAndEnergy(rij, Rij, r0, fc);

        // Total.
        auto energy{pairSR.first + pairHP.first};
        auto fi{pairSR.second + pairHP.second};

        // Done.
        logger.trace("Leaving");
        return std::move(std::make_pair(energy, fi));
    }
}