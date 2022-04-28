/*
 * Author: André H. Juffer.
 * Created on 21/12/2021, 14:30.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef UTIL_DIRECTION_HPP
#define UTIL_DIRECTION_HPP

namespace simploce {

    /**
     * Specifies a direction along which something is applied.
     * x : Applied along the x-direction.
     * y : Applied along the y-direction.
     * z : Applied along the z-direction.
     */

    class Direction {
    public:

        /**
         * Applied along the x-direction.
         */
        static Direction X;

        /**
         * Applied along the y-direction.
         */
        static Direction Y;

        /**
         * Applied along the y-direction.
         */
        static Direction Z;

        /**
         * Returns direction.
         * @param value Direction value. One of 'x', 'y', or 'z'.
         * @return Direction.
         */
        static Direction valueOf(char value = 'z');

        /**
         * Equality operator.
         * @param direction Direction.
         * @return Result
         */
        bool operator == (const Direction& direction) const;

    private:

        explicit Direction(char value);

        char value_;

    };

}

#endif //UTIL_DIRECTION_HPP
