/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:55 PM
 */

#include "simploce/simulation/no-bc.hpp"
#include "simploce/conf/s-conf.hpp"

namespace simploce {

    dist_vect_t 
    NoBoundaryCondition::apply(const position_t& ri,
                               const position_t& rj) const {
        //return std::move(dist_vect_t{ri - rj});
        return boundary_condition_impl::apply(ri, rj);
    }
    
    position_t 
    NoBoundaryCondition::placeInside(const position_t& r_out) const {
        //return r_out;
        return boundary_condition_impl::placeInside(r_out);
    }

    position_t
    NoBoundaryCondition::apply(const position_t &r) const {
        return boundary_condition_impl::apply(r);
    }

    velocity_t
    NoBoundaryCondition::apply(const simploce::velocity_t &v,
                               const position_t& r) const {
        return boundary_condition_impl::apply(v, r);
    }

    void
    NoBoundaryCondition::applyToVelocities(const simploce::pg_ptr_t &particleGroup) const {
        boundary_condition_impl::applyToVelocities(particleGroup);
    }

}
