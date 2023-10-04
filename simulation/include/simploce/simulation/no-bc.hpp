/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 5:52 PM
 */

#ifndef NO_BC_HPP
#define NO_BC_HPP

#include "bc.hpp"

namespace simploce {

    /**
     * No boundary conditions are applied.
     */
    struct NoBoundaryCondition : public boundary_condition {

        dist_vect_t apply(const position_t& ri, const position_t& rj) const override;
        
        position_t placeInside(const position_t& r_out) const override;

        /**
         * No effect.
         */
        velocity_t apply(const simploce::velocity_t &v,
                         const position_t& r) const override;

        /**
         * No effect.
         */
        void applyToVelocities(const pg_ptr_t& particleGroup) const override;

    };
}



#endif /* NO_BC_HPP */

