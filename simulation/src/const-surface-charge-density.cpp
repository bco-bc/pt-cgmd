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
#include <cmath>

namespace simploce {

    ConstantSurfaceChargeDensity::ConstantSurfaceChargeDensity(srf_charge_density_t sigma,
                                                               FlatSurface flatSurface,
                                                               real_t eps_r,
                                                               bc_ptr_t bc) :
        sigma_{sigma}, flatSurface_{flatSurface}, eps_r_{eps_r}, bc_{bc} {
    }

    std::pair<energy_t, force_t>
    ConstantSurfaceChargeDensity::operator () (const p_ptr_t& particle) {
        auto r = particle->position();
        return std::move(ConstantSurfaceChargeDensity::forceAndEnergy(sigma_,
                                                                      flatSurface_,
                                                                      eps_r_,
                                                                      bc_,
                                                                      particle->position(),
                                                                      particle->charge()));
    }

    std::pair<energy_t, force_t>
    ConstantSurfaceChargeDensity::forceAndEnergy(srf_charge_density_t sigma,
                                                 const FlatSurface& flatSurface,
                                                 real_t eps_r,
                                                 const bc_ptr_t& bc,
                                                 const position_t& r,
                                                 const charge_t& q) {
        auto pair = flatSurface.distanceTo(r);
        auto R = pair.first;
        energy_t energy{-sigma() * R() * q()/ (2.0 * units::mu<real_t>::E0 * eps_r)};
        real_t dUrdR =  -sigma() / (2.0 * units::mu<real_t>::E0 * eps_r);
        dist_vect_t unitVector = flatSurface.unitVectorPerpendicularTo();
        force_t f{};
        for (int k = 0; k != 3; ++k) {
            f[k] = -dUrdR * unitVector[k] * q();
        }
        return std::move(std::make_pair(energy, f));
    }
}