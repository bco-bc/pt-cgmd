/*
 * File:   math-constants.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 3:48 PM
 */

#ifndef MATH_CONSTANTS_HPP
#define MATH_CONSTANTS_HPP

#include <cmath>
#include <limits>

namespace simploce {
    namespace math {

        /**
         * Holds various values for some mathematical constants.
         * @param V Value type (e.g. double).
         */
        template<typename V>
        struct constants {

            /**
             * Value for PI.
             * @see <a href="http://en.wikipedia.org/wiki/Pi">PI and Wikipedia</a>
             */
            static const V PI;

            /**
             * Value for the number E (or e).
             * @see <a href="http://en.wikipedia.org/wiki/E_%28mathematical_constant%29">E at Wikipedia</a>
             */
            static const V E;

            /**
             * Value for the Euler constant.
             * @see <a href="http://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant">
             *   Euler constant at Wikipedia
             * </a>
             */
            static const V GAMMA;

            /**
             * A very large value.
             */
            static const V VERY_LARGE;

            /**
             * A large value
             */
            static const V LARGE;

            /**
             * A very small value.
             */
            static const V VERY_SMALL;

            /**
             * A small value.
             */
            static const V SMALL;

        };

        template<typename V>
        const V constants<V>::PI = std::acos(-1.L);

        template<typename V>
        const V constants<V>::E = 2.71828182845904523536028747135266249775724709369995;

        template<typename V>
        const V constants<V>::GAMMA = 0.5772156649015328606065120900824024310421;

        template<typename V>
        const V constants<V>::VERY_LARGE = std::numeric_limits<V>::max();

        template<typename V>
        const V constants<V>::LARGE = 1.0e+16;

        template<typename V>
        const V constants<V>::VERY_SMALL = std::numeric_limits<V>::epsilon();

        template<typename V>
        const V constants<V>::SMALL = 1.0e-16;

    }
}

#endif /* MATH_CONSTANTS_HPP */

