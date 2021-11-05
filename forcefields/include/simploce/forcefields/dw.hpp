/*
 * File: dw.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 16, 2021., 1:00 PM.
 */

#ifndef DW_HPP
#define DW_HPP

#include "ext-potential.hpp"

namespace simploce {
    /**
     * Interface for double well potential, U(x) where x is the location in the double
     * well.
     */
    template <typename P>
    struct dw : public ext_potential<P> {

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = typename ext_potential<P>::particle_ptr_t;

        /**
         * Value of U(x)
         * @param x Location.
         * @return U(x)
         */
        virtual real_t U(real_t x) const = 0;

        /**
         * First derivative of U(x) at x.
         * @param x Location.
         * @return First derivative.
         */
        virtual real_t dUdx(real_t x) const = 0;

        /**
         * Second derivative of U(x) at x.
         * @param x Location.
         * @return Second derivative.
         */
        virtual real_t dU2dx2(real_t x) const = 0;
    };
}

#endif //DW_HPP
