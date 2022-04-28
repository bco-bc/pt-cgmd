/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 13:37.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef UTIL_PLANE_HPP
#define UTIL_PLANE_HPP

#include <string>

namespace simploce {

    /**
     * One of the three planes (xy, yx, and xz plane) in a cartesian coordinate system.
     */
    class Plane {
    public:

        /**
         * XY-plane
         */
        static Plane XY;

        /**
         * YX-plane
         */
        static Plane YZ;

        /**
         * XZ-plane.
         */
        static Plane ZX;

        /**
         * Returns plane.
         * @param value Plane value. One of "xy", "yz", and "xz". Any other value returns Plane::XY.
         * @return Plane.
         */
        static Plane valueOf(const std::string& value = "xy");

        /**
         * Equality operator.
         * @param plane Plane.
         * @return Result.
         */
        bool operator == (const Plane& plane);

        /**
         * Returns string representation.
         * @return String representation.
         */
        std::string toString() const;

    private:

        explicit Plane(std::string value);

        std::string value_;
    };

}

#endif //UTIL_PLANE_HPP
