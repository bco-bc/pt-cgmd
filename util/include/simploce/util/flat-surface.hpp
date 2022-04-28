/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 13:28.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef UTIL_FLAT_SURFACE_HPP
#define UTIL_FLAT_SURFACE_HPP

#include "../types/u-types.hpp"
#include "plane.hpp"
#include <iostream>

namespace simploce {

    /**
     * Parallel to a given plane and located at the distance away from that plane. Coordinates are
     * in the positive quadrant of the cartesian coordinate system.
     */
    class FlatSurface {
    public:

        /**
         * Constructor. The default values places the surface parallel to the xy-plane at z = 0.
         * @param plane Plane
         * @param distanceToPlane Distance to plane.
         */
        explicit FlatSurface(Plane plane = Plane::XY, const dist_t& distanceToPlane = dist_t{0.0});

        /**
         * Returns the distance between a given point and the surface.
         * @param r Point.
         * @return Distance.
         */
        std::pair<dist_t, dist_vect_t> distanceTo(const position_t& r) const;

        /**
         * Returns an unit vector perpendicular to the surface. The unit normal is pointing towards
         * the positive x/y/z axis.
         * @return Unit vector.
         */
        dist_vect_t unitVectorPerpendicularTo() const;

        /**
         * Returns string representation.
         * @return String.
         */
        std::string toString() const;

    private:

        friend std::ostream& operator <<(std::ostream& stream, const FlatSurface& flatSurface);

        Plane plane_;
        dist_t distanceToPLane_;
        dist_vect_t unitVector_;

        size_t coordinate_ = 2;
    };

    /**
     * Writes flat surface to an output stream.
     * @param stream Output stream.
     * @param flatSurface Flat surface.
     * @return Output stream.
     */
    std::ostream& operator <<(std::ostream& stream, const FlatSurface& flatSurface);

}

#endif //UTIL_FLAT_SURFACE_HPP
