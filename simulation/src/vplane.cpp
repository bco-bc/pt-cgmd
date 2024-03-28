/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 28 February 2024. Adapted from original version (1996).
 * @see https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(199612)17:16%3C1783::AID-JCC1%3E3.0.CO;2-J
 */

#include "simploce/potentials/vplane.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/math-constants.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <cmath>
#include <stdexcept>

namespace simploce {
    namespace virtual_plane {

        /**
         * Calculates the interaction energy of the given charge with the given virtual plane.
         * Energy only, no forces as of yet.
         * @param planeLocation Distance of the plane to the xy-plane at z = 0.
         * @param sigma Plane's surface charge density.
         * @param Lx Length of box in the x-direction.
         * @param eps_r Relative permittivity.
         * @param r Particle position -inside- the box.
         * @param Q Particle charge.
         * @return Energy and forces.
         */
        static std::pair<energy_t, force_t>
        interaction(const dist_t &planeLocation,
                    const srf_charge_density_t &sigma,
                    const length_t &Lx,
                    real_t eps_r,
                    const position_t& r,
                    const charge_t &Q) {
            static util::Logger logger{"simploce::virtual_plane::interaction()"};
            logger.trace("Entering.");

            static auto PI = math::constants<real_t>::PI;
            static auto E0 = units::mu<real_t>::E0;
            static auto PI_eps_r_E0 = PI * eps_r * E0;
            static auto two_eps_r_E0 = 2.0 * eps_r * E0;
            static auto constant_at_0 =
                    Lx() * std::log(std::tan(3.0 * PI / 8.0)) / PI_eps_r_E0;  // At R = 0.0.
            static auto two_over_PI_eps_r_E0 = 2.0 / PI_eps_r_E0;
            static auto Lx2 = Lx() * Lx();
            static auto quarter_Lx2 = Lx2 / 4.0;

            // Calculate integral from 0 to 1/4*PI.
            static real_t a = 0;
            static real_t b = PI / 4.0;
            auto R = std::fabs(planeLocation() - r[2]);
            auto R2 = R * R;
            auto integrand = [R2] (real_t x) {
                auto cos_x = std::cos(x);
                auto cos_2_x = cos_x * cos_x;
                auto f_x = std::sqrt((quarter_Lx2 + R2 * cos_2_x) / cos_2_x);
                return f_x;
            };
            auto integral = util::integrate(integrand, a, b);

            // Energy.
            auto energy = -Q() * sigma() * (two_over_PI_eps_r_E0 * integral - constant_at_0);

            logger.trace("Leaving.");
            return std::move(std::make_pair(energy_t{energy}, force_t{}));
        }

    }

    VirtualPlane::VirtualPlane(box_ptr_t box,
                               bc_ptr_t bc,
                               dist_t location,
                               real_t eps_r) :
            box_{std::move(box)}, bc_{std::move(bc)}, location_{location}, eps_r_{eps_r},
            sigma_{0.0} {
        if (!box_) {
            throw std::domain_error("Box must be provided.");
        }
        if (!bc_) {
            throw std::domain_error("Boundary conditions must be provided.");
        }
        if (eps_r_ <= 0) {
            throw std::domain_error("The relative permittivity must be a positive number.");
        }
    }

    srf_charge_density_t
    VirtualPlane::surfaceChargeDensity() const {
        return sigma_;
    }

    dist_t
    VirtualPlane::location() const {
        return location_;
    }

    std::pair<energy_t, force_t>
    VirtualPlane::operator()(const simploce::p_ptr_t &particle) const {
        static auto Lx = box_->lengthX();

        // Particle charge and position.
        auto Q = particle->charge();
        auto r = particle->position();
        r = bc_->placeInside(r);

        // Interaction with this plane.
        auto result = virtual_plane::interaction(location_,
                                                 sigma_,
                                                 Lx,
                                                 eps_r_,
                                                 r,
                                                 Q);

        // Done.
        return std::move(result);
    }

    void
    VirtualPlane::reset(srf_charge_density_t sigma) {
        sigma_ = sigma;
    }

    std::ostream&
    operator << (std::ostream& stream, const VirtualPlane& plane) {
        stream << plane.location() << plane.surfaceChargeDensity();
        return stream;
    }
}