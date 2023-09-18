/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 12:25.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/voltage.hpp"
#include "simploce/potentials/const-surface-charge-density.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include <memory>

namespace simploce {

    /*
     * This implementation creates a flat surface with a uniform surface density located at one of the faces of the
     * simulation box. The electric potential and field is obtained by placing a particle next to it.
     * The electric potential's reference point is at surface (not at infinity).
     */

    static std::shared_ptr<ConstantSurfaceChargeDensity> chargedSurface_{};

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

    static std::shared_ptr<ConstantSurfaceChargeDensity>
    createChargedSurface(const Direction& direction,
                         voltage_t voltage,
                         dist_t distance,
                         real_t eps_r,
                         const bc_ptr_t& bc,
                         bool mesoscopic) {
        FlatSurface flatSurface = create_(direction);
        auto e0 = units::mu<real_t>::E0;
        real_t epd = units::mu<real_t>::V_to_kJ_mol_e * voltage();           // Conversion V to kJ/(mol e)
        srf_charge_density_t sigma = -epd * 2.0 * e0 * eps_r / distance();   // In e/nm²
        return std::make_shared<ConstantSurfaceChargeDensity>(sigma, flatSurface, eps_r, bc, mesoscopic);
    }

    Voltage::Voltage(voltage_t voltage,
                     dist_t distance,
                     real_t eps_r,
                     const bc_ptr_t& bc,
                     const Direction& direction,
                     bool mesoscopic) {
        chargedSurface_ = createChargedSurface(direction, voltage, distance, eps_r, bc, mesoscopic);
    }

    std::pair<energy_t, force_t>
    Voltage::operator () (const p_ptr_t& particle) {
        static util::Logger logger{"simploce::Voltage::operator () "};
        logger.trace("Entering");

        logger.debug("This external potential is included.");

        auto result = chargedSurface_->forceAndEnergy(particle->position(), particle->charge());

        logger.trace("Leaving.");
        return std::move(result);
    }

}

