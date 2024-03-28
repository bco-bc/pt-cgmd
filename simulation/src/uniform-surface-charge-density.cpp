/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#include "simploce/potentials/uniform-surface-charge-density.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/conf/s-conf.hpp"
#include <utility>

namespace simploce {

    UniformSurfaceChargeDensity::UniformSurfaceChargeDensity(srf_charge_density_t sigma,
                                                             FlatSurface flatSurface,
                                                             real_t eps_r,
                                                             bc_ptr_t bc,
                                                             dist_t delta,
                                                             bool mesoscopic) :
            sigma_{sigma}, flatSurface_{std::move(flatSurface)}, eps_r_{eps_r}, bc_{std::move(bc)},
            delta_{delta}, mesoscopic_{mesoscopic} {
        util::Logger logger{"simploce::UniformSurfaceChargeDensity::UniformSurfaceChargeDensity()"};
        logger.trace("Entering");

        logger.info(std::to_string(sigma_()) + ": Surface charge density.");
        logger.info(flatSurface_.toString() + ": Flat surface description.");
        logger.info(std::to_string(eps_r_) + ": Relative permittivity.");
        logger.info(std::to_string(delta_()) + ": Width Stern layer.");
        logger.info(std::to_string(mesoscopic_) + ": Mesoscale?");

        logger.trace("Leaving.");
    }

    std::pair<energy_t, force_t>
    UniformSurfaceChargeDensity::operator () (const p_ptr_t& particle) const {
        static util::Logger logger("simploce::UniformSurfaceChargeDensity::operator () ()");
        logger.trace("Entering");

        logger.debug("This external potential is included.");
        auto r = particle->position();
        auto Q = particle->charge();
        auto radius = particle->spec()->radius();
        auto result = this->forceAndEnergy(r, radius, Q);

        logger.trace("Leaving");
        return std::move(result);
    }

    std::pair<energy_t, force_t>
    UniformSurfaceChargeDensity::forceAndEnergy(const position_t& r,
                                                const radius_t& radius,
                                                const charge_t& Q) const {
        static util::Logger logger("simploce::UniformSurfaceChargeDensity::forceAndEnergy");
        logger.trace("Entering.");

        auto ri = bc_->placeInside(r);
        auto pair = flatSurface_.distanceTo(ri);
        auto R = pair.first;
        logger.debug(std::to_string(R()) + ": Distance to surface.");
        logger.debug(util::to_string(Q) + ": Charge value.");
        if ( R <= (delta_ + radius) ) {
            return std::move(std::make_pair(conf::LARGE, force_t{0.0, 0.0, 0.0}));
        }
        auto E0 = mesoscopic_ ? 1.0 : units::mu<real_t>::E0;
        energy_t energy{-sigma_() * R() * Q() / (2.0 * E0 * eps_r_)};
        real_t dUrdR =  -sigma_() / (2.0 * E0 * eps_r_);
        dist_vect_t unitVector = flatSurface_.unitVectorPerpendicularTo();
        force_t f{};
        for (int k = 0; k != 3; ++k) {
            f[k] = -dUrdR * unitVector[k] * Q();
        }
        logger.debug(std::to_string(energy()) + ": Energy.");

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, f));
    }

}