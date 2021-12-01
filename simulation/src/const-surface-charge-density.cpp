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

    /**
     * Returns distance and distance vector to the surface.
     * @param plane Location of the surface.
     * @param r Position of some particle.
     * @param bc Boundary condition.
     * @return Distance (length of the distance vector) and distance vector.
     */
    static std::pair<dist_t, dist_vect_t>
    distanceToSurface(ConstantSurfaceChargeDensity::PLANE plane,
                      const position_t& r,
                      const bc_ptr_t& bc) {
        dist_vect_t r_in = bc->placeInside(r);
        if (plane == ConstantSurfaceChargeDensity::xy ) {
            return std::move(std::make_pair(std::fabs(r_in[2]), dist_vect_t{0.0, 0.0, r_in[2]}));
        } else if (plane == ConstantSurfaceChargeDensity::yz) {
            return std::move(std::make_pair(std::fabs(r_in[0]), dist_vect_t{r_in[0], 0.0, 0.0}));
        } else {
            return std::move(std::make_pair(std::fabs(r_in[1]), dist_vect_t{0.0, r_in[1], 0.0}));
        }
    }

    ConstantSurfaceChargeDensity::PLANE
    ConstantSurfaceChargeDensity::valueOf(std::string value) {
        if ( value == "xy") {
            return ConstantSurfaceChargeDensity::xy;
        } else if (value == "yz") {
            return ConstantSurfaceChargeDensity::yz;
        } else {
            return ConstantSurfaceChargeDensity::zx;
        }
    }

    ConstantSurfaceChargeDensity::ConstantSurfaceChargeDensity(srf_charge_density_t sigma,
                                                               real_t eps_r,
                                                               PLANE plane,
                                                               bc_ptr_t bc) :
        sigma_{sigma}, eps_r_{eps_r}, plane_(plane), bc_{std::move(bc)} {
    }

    std::pair<energy_t, force_t>
    ConstantSurfaceChargeDensity::operator () (const p_ptr_t& particle) {
        auto r = particle->position();
        return std::move(ConstantSurfaceChargeDensity::forceAndEnergy(sigma_,
                                                                      plane_,
                                                                      eps_r_,
                                                                      bc_,
                                                                      particle->position(),
                                                                      particle->charge()));
    }

    std::pair<energy_t, force_t>
    ConstantSurfaceChargeDensity::forceAndEnergy(srf_charge_density_t sigma,
                                                 PLANE plane,
                                                 real_t eps_r,
                                                 const bc_ptr_t& bc,
                                                 const position_t& r,
                                                 const charge_t& q) {
        auto pair = distanceToSurface(plane, r, bc);
        energy_t energy{-sigma() * pair.first() / (2.0 * units::mu<real_t>::E0 * eps_r)};
        real_t dUrdR =  -sigma() / (2.0 * units::mu<real_t>::E0 * eps_r);
        dist_vect_t unitVector = pair.second / pair.first;
        force_t f{};
        for (int k = 0; k != 3; ++k) {
            f[k] = -dUrdR * unitVector[k] * q();
        }
        return std::move(std::make_pair(energy, f));
    }
}