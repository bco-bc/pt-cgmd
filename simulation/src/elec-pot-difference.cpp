/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 12:25.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/elec-pot-difference.hpp"
#include "simploce/potentials/const-surface-charge-density.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include <memory>

namespace simploce {

    static ConstantSurfaceChargeDensity::PLANE plane_{};
    static srf_charge_density_t sigma_{};                // e/nm^2
    static real_t eps_r_{};

    static ConstantSurfaceChargeDensity::PLANE
    surfaceLocation(ElectricPotentialDifference::DIRECTION direction) {
        if ( direction == ElectricPotentialDifference::x) {
            return ConstantSurfaceChargeDensity::yz;
        } else if (direction == ElectricPotentialDifference::y) {
            return ConstantSurfaceChargeDensity::zx;
        } else {
            return ConstantSurfaceChargeDensity::xy;
        }
    }

    ElectricPotentialDifference::DIRECTION
    ElectricPotentialDifference::valueOf(char value) {
        if (value == 'x') {
            return ElectricPotentialDifference::x;
        } else if (value == 'y') {
            return ElectricPotentialDifference::y;
        } else {
            return ElectricPotentialDifference::z;
        }
    }

    ElectricPotentialDifference::ElectricPotentialDifference(el_pot_diff deltaV,
                                                             dist_t distance,
                                                             real_t eps_r,
                                                             bc_ptr_t bc,
                                                             DIRECTION direction) :
        bc_{std::move(bc)} {
        plane_ = surfaceLocation(direction);
        eps_r_ = eps_r;
        auto e0 = units::mu<real_t>::E0;
        real_t epd = units::mu<real_t>::V_to_kJ_mol_e * deltaV();  // Conversion V to kJ/(mol e)
        sigma_ = -epd * 2.0 * e0 * eps_r_ / distance();
    }

    std::pair<energy_t, force_t>
    ElectricPotentialDifference::operator () (const p_ptr_t& particle) {
        return std::move(ConstantSurfaceChargeDensity::forceAndEnergy(sigma_,
                                                                         plane_,
                                                                         eps_r_,
                                                                         bc_,
                                                                         particle->position(),
                                                                         particle->charge()));
    }

}

