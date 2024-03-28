/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on February 21, 2024.
 */

#include "simploce/simulation/bc-impl.hpp"

namespace simploce {

    boundary_condition_impl::boundary_condition_impl() = default;

    boundary_condition_impl::~boundary_condition_impl() noexcept = default;

    dist_vect_t
    boundary_condition_impl::apply(const position_t &ri,
                                   const position_t &rj) const {
        return dist_vect_t{ri - rj};
    }

    velocity_t
    boundary_condition_impl::apply(const velocity_t &v,
                                   const position_t &r) const {
        return v;
    }

    position_t
    boundary_condition_impl::placeInside(const position_t &r_out) const {
        return r_out;
    }

    position_t
    boundary_condition_impl::apply(const position_t &r) const {
        return r;
    }

    void
    boundary_condition_impl::applyToVelocities(const pg_ptr_t &particleGroup) const {
    }
}
