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
#include "simploce/util/logger.hpp"
#include <memory>

namespace simploce {

    FlatSurface
    create_(const Direction& direction) {
        if (direction == Direction::X) {
            return FlatSurface{Plane::YZ};
        } else if (direction == Direction::Y) {
            return FlatSurface{Plane::ZX};
        } else {
            return FlatSurface{Plane::XY};
        }
    }

    static FlatSurface flatSurface_{};
    static srf_charge_density_t sigma_{};                // e/nm^2

    ElectricPotentialDifference::ElectricPotentialDifference(el_pot_diff deltaV,
                                                             dist_t distance,
                                                             real_t eps_r,
                                                             bc_ptr_t bc,
                                                             const Direction& direction) :
        bc_{std::move(bc)}, eps_r_{eps_r} {
        flatSurface_ = create_(direction);
        auto e0 = units::mu<real_t>::E0;
        real_t epd = units::mu<real_t>::V_to_kJ_mol_e * deltaV();  // Conversion V to kJ/(mol e)
        sigma_ = -epd * 2.0 * e0 * eps_r_ / distance();
    }

    std::pair<energy_t, force_t>
    ElectricPotentialDifference::operator () (const p_ptr_t& particle) {
        static util::Logger logger{"simploce::ElectricPotentialDifference::operator () "};
        logger.debug("This external potential is included.");
        return std::move(ConstantSurfaceChargeDensity::forceAndEnergy(sigma_,
                                                                      flatSurface_,
                                                                      eps_r_,
                                                                      bc_,
                                                                      particle->position(),
                                                                      particle->charge()));
    }

}

