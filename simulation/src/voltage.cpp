/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 12:25.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/voltage.hpp"
#include "simploce/potentials/uniform-surface-charge-density.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <memory>

namespace simploce {

    /*
     * This implementation creates a flat surface with a uniform surface density. The surface is located at one of
     * the faces of the simulation box depending on the direction of the electric field. For instance, if the
     * electric field should point in the positive direction, then the charged surface is placed at the xy-plane
     * at z=0. The reference point for the electric potential is at z=0 (and not at infinity).
     */

    static std::shared_ptr<UniformSurfaceChargeDensity> chargedSurface_{};

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

    static std::shared_ptr<UniformSurfaceChargeDensity>
    createChargedSurface(const Direction& direction,
                         voltage_t voltage,
                         dist_t distance,
                         real_t eps_r,
                         const bc_ptr_t& bc,
                         bool mesoscopic) {
        FlatSurface flatSurface = create_(direction);
        auto e0 = mesoscopic ? 1.0 : units::mu<real_t>::E0;
        auto V_to_kJ_mol_e = mesoscopic ? 1.0 : units::mu<real_t>::V_to_kJ_mol_e;
        real_t epd = V_to_kJ_mol_e * voltage();
        srf_charge_density_t sigma = -epd * 2.0 * e0 * eps_r / distance();
        return std::make_shared<UniformSurfaceChargeDensity>(sigma, flatSurface, eps_r, bc, mesoscopic);
    }

    Voltage::Voltage(voltage_t voltage,
                     dist_t distance,
                     bc_ptr_t bc,
                     real_t eps_r,
                     bool mesoscopic) :
        e0_{el_field_t {0.0, 0.0, -voltage() / distance()}}, bc_{std::move(bc)}, eps_r_{eps_r}, mesoscopic_{mesoscopic} {
    }

    Voltage::Voltage(simploce::el_field_t e0, simploce::bc_ptr_t bc, real_t eps_r, bool mesoscopic) :
        e0_{e0}, bc_{std::move(bc)}, eps_r_{eps_r}, mesoscopic_{mesoscopic} {
        util::Logger logger{"simploce::Voltage::Voltage()"};
        logger.trace("Entering.");

        if ( !bc_ ) {
             throw std::domain_error("Boundary conditions must be provided.");
        }
        if (eps_r_ <= 0) {
            throw std::domain_error("Relative permittivity must be a positive number.");
        }

        logger.info(util::to_string(e0_) + ": External static homogeneous electric field.");
        logger.info(util::to_string(eps_r_) + ": Relative permittivity.");
        logger.info(std::to_string(mesoscopic_) + ": Mesoscopic simulation?");

        logger.trace("Leaving.");
    }

    std::pair<energy_t, force_t>
    Voltage::operator () (const p_ptr_t& particle) const {
        static util::Logger logger{"simploce::Voltage::operator()"};
        logger.trace("Entering.");

        logger.debug("This external potential is included.");
        auto r = bc_->placeInside(particle->position());
        r = particle->position();
        auto Q = particle->charge();
        energy_t energy{-Q() * inner<real_t>(r, e0_) / eps_r_};
        force_t f{0.0, 0.0, 0.0};
        for (int k = 0; k != 3; ++k) {
            f[k] = Q() * e0_[k] / eps_r_;
        }

        logger.trace("Leaving.");
        return std::move(std::make_pair(energy, f));
    }

}

