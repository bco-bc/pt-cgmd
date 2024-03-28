/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on February 21, 2024.
 */

#ifndef SIMULATION_BC_IMPL_HPP
#define SIMULATION_BC_IMPL_HPP

#include "bc.hpp"

namespace simploce {

    /**
     * Default implementation from which boundary conditions may be derived.
     */
    struct boundary_condition_impl : public boundary_condition {

        boundary_condition_impl();

        ~boundary_condition_impl() noexcept override;

        /**
         * @return ri - rj
         */
        dist_vect_t
        apply(const position_t& ri,
              const position_t& rj) const override;

        /**
         * No effect.
         * @return r_out
         */
        position_t
        placeInside(const position_t& r_out) const override;

        /**
         * No effect.
         * @return r
         */
        position_t
        apply(const position_t& r) const override;


        /**
         * No effect
         * @return v
         */
        velocity_t
        apply(const velocity_t &v,
              const position_t& r) const override;

        /**
         * No effect.
         */
        void
        applyToVelocities(const pg_ptr_t& particleGroup) const override;

    };
}

#endif //SIMULATION_BC_IMPL_HPP
