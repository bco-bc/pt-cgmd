/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#include "simploce/potentials/const-surface-charge-density.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include <cmath>
#include <utility>

namespace simploce {

    ConstantSurfaceChargeDensity::ConstantSurfaceChargeDensity(srf_charge_density_t sigma,
                                                               FlatSurface flatSurface,
                                                               real_t eps_r,
                                                               bc_ptr_t bc,
                                                               bool mesoscopic) :
        sigma_{sigma}, flatSurface_{std::move(flatSurface)}, eps_r_{eps_r}, bc_{std::move(bc)},
        mesoscopic_{mesoscopic} {
    }

    std::pair<energy_t, force_t>
    ConstantSurfaceChargeDensity::operator () (const p_ptr_t& particle) {
        static util::Logger logger("simploce::ConstantSurfaceChargeDensity::operator () ()");
        logger.trace("Entering");

        auto r = particle->position();
        auto Q = particle->charge();
        auto result = this->forceAndEnergy(r, Q);

        logger.trace("Leaving");
        return std::move(result);
    }

    std::pair<energy_t, force_t>
    ConstantSurfaceChargeDensity::forceAndEnergy(const position_t& ro, const charge_t& q, real_t df) {
        auto r = bc_->placeInside(ro);
        auto pair = flatSurface_.distanceTo(r);
        auto R = pair.first;
        energy_t energy{-sigma_() * R() * q()/ (2.0 * units::mu<real_t>::E0 * eps_r_)};
        real_t dUrdR =  -sigma_() / (2.0 * units::mu<real_t>::E0 * eps_r_);
        dist_vect_t unitVector = flatSurface_.unitVectorPerpendicularTo();
        force_t f{};
        for (int k = 0; k != 3; ++k) {
            f[k] = -dUrdR * unitVector[k] * q();
        }
        return std::move(std::make_pair(energy, f));
    }
}