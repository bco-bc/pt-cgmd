/*
 * File: dw.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 26, 2021., 11:00 PM.
 */

#ifndef FORCEFIELDS_STATIC_UNIFORM_HPP
#define FORCEFIELDS_STATIC_UNIFORM_HPP

#include "ext-potential.hpp"
#include "simploce/types/u-types.hpp"
#include "simploce/units/units-mu.hpp"
#include <cmath>
#include <memory>

namespace simploce {

    /**
     * The electric potential and field due to an infinite surface parallel to yz-plane
     * at x=0 with a static uniform surface charge density. The reference point for the
     * electric potential is at x=0.
     * @tparam P
     */
    template <typename P>
    class StaticUniform : public ext_potential<P> {
    public:

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = typename ext_potential<P>::particle_ptr_t;

        /**
         * Constructor
         * @param sigma Surface charge density.
         * @param e Relative permittivity.
         */
        StaticUniform(const srf_charge_density_t& sigma,
                      const rel_perm_t& eps_r) : ext_potential<P>{}, sigma_{sigma}, eps_r_{eps_r} {};

        /**
         * Energy.
         * @param p Particle.
         * @return energy in kJ/mol.
         */
        energy_t energy(const particle_ptr_t& p) const override;

        /**
         * Adds force due to the surface charge density.
         * @param p Particle.
         * @return Energy in kJ/mol.
         */
        energy_t force(particle_ptr_t& p) const override;

    private:

        real_t U_(real_t x) const {
            return -sigma_() * x / (2.0 * units::mu<real_t>::E0 * eps_r_());
        }

        real_t dUdx_() const {
            return -sigma_() / (2.0 * units::mu<real_t>::E0 * eps_r_());
        }

        srf_charge_density_t sigma_;
        rel_perm_t eps_r_;

    };

    template <typename P>
    energy_t
    StaticUniform<P>::energy(const particle_ptr_t& p) const {
        position_t r = p->position();
        distance_t R = std::fabs(r[0]);
        return p->charge()() * U_(R());
    }


    template <typename P>
    energy_t
    StaticUniform<P>::force(particle_ptr_t& p) const {
        position_t r = p->position();
        auto x = r[0];
        distance_t R = std::fabs(x);
        auto dudx = dUdx_();
        auto f_x = -p->charge()() * dudx * x / std::fabs(x);
        energy_t energy{R() * p->charge()() * dudx};
        force_t f{f_x, 0.0, 0.0};
        f += p->force();
        p->force(f);
        return energy;
    }

}

#endif //FORCEFIELDS_STATIC_UNIFORM_HPP
