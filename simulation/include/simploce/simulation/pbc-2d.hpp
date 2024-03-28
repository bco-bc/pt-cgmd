/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on February 1, 2024.
 */


#ifndef SIMULATION_PBC_2D_HPP
#define SIMULATION_PBC_2D_HPP

#include "bc-impl.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * Periodic boundary conditions with the nearest image approximation applied to two directions,
     * e.g. x and y. The remaining third coordinate is applied to a given position such that the
     * coordinate (e.g., z) always falls in the range [0, box_k], where box_k is the box side in that
     * direction (e.g., box size in the z-direction), i.e. it is "reinserted".
     */
    class PBC_2D : public boundary_condition_impl {
    public:

        /**
         * Constructor.
         * @param box Box
         * @param d1 First direction.
         * @param d2 Second direction.
         * @param reinsert Coordinate that should be "reinserted".
         */
        explicit PBC_2D(box_ptr_t box,
                        Direction d1 = Direction::X,
                        Direction d2 = Direction::Y,
                        Direction reinsert = Direction::Z);

        dist_vect_t
        apply(const position_t& ri,
              const position_t& rj) const override;

        position_t
        placeInside(const position_t& r_out) const override;

        /**
         * No effect.
         */
        velocity_t apply(const simploce::velocity_t &v,
                         const position_t& r) const override;

        /**
         * Randomly determineStateChanges the "keepInBox" coordinate.
         * @param r Position.
         * @return Reinserted position.
         */
        position_t
        apply(const position_t& r) const override;

        /**
         * No effect.
         */
        void applyToVelocities(const pg_ptr_t& particleGroup) const override;

    private:

        box_ptr_t box_;
        int directions_[2];
        int reinsert_;
    };
}

#endif //SIMULATION_PBC_2D_HPP
