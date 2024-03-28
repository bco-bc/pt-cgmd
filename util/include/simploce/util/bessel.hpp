/*
 * Author: André H. Juffer.
 * Created on 5/2/2024
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef UTIL_BESSEL_HPP
#define UTIL_BESSEL_HPP

#include "../types/u-types.hpp"

namespace simploce {
    namespace math {

        /*
         * Some implementations of methods to compute Bessel function of the first, second,
         * third kind, and Bessel functions of imaginary argument.
         * See for an introduction: N.N. Lebedev, "Special functions and their
         * applications", Dover Publications, Inc., New York, 1972.
         *
         * Real numbers are assumed throughout.
         */

        /*
         * This method decides on the basis on the value of the function argument x which
         * method should be used for computing the Bessel function J0(x) of the first kind
         * of order 0.
         */
        real_t
        Bessel_J0(real_t x);

        /* Bessel function Jn(z) of the first kind of order n, where n is of nonnegative
         * integral order. This is a straightforward implementation of the series for
         * Jn(z) and it will converges in the whole complex plane.
         * It will only work accurately for small -8<=z<=8 and small n (0,+1,-1).
         * For other cases, one must use the recurrence relations (see below). But this
         * function may provide J0 and J1 to start the recurrence.
         */
        real_t
        B1_Jn_S(size_t n,
                real_t &x);    // S stands for Series.

        /*
         * Bessel function J0(x) determined by several polynomial approximations (Handbook
         * of Mathematical functions, Abramowitz and Stengun (AS), Dover Publications, Inc.
         * New York, 1972)
         */

        // The next approximation is valid for -3<=x<=3 with an |error| < 5e-08.
        real_t
        B1_J0_PA_AS941(real_t x);    // PA stands for polynomial approximation.
                                     // 941 stands for formula 9.4.1 in AS,
                                     // page 369.

        // The next approximation is valid for x>=3 with an |error| < ...
        real_t
        B1_J0_PA_AS943(real_t x);    // AS, 9.4.3, page 369.

        /*
         * Bessel function Jn(z) of the first kind of order n, where n is of nonnegative
         * integral order. This makes use of the recurrence relation:
         * z*J_{n-1}(z)+z*J_{n+1}(z)=2*n*Jn(z). Initial values for J0 and J1 are
         * computed by other functions.
         *
        real_t
        B1_Jn_RR(size_t n, real_t x);     // RR stands for recurrence relation.
        */

        /*
         * Modified Bessel function of the third kind In(z) and MacDonald's function
         * Kn(z). These functions are also called Bessel functions of imaginary
         * argument. All functions use a polynomial approximation.
         */

        real_t
        Bessel_K0(real_t x);          // Valid for all x>0.

        real_t
        Bessel_K1(real_t x);          // Valid for all x>0.

        real_t
        Bessel_I0(real_t x);          // Valid for x>-3.75.

        real_t
        Bessel_I1(real_t x);          // Valid for x>-3.75.

        real_t
        B3M_I0_PA_AS981(real_t x);    // |error| < 1.6e-07 for -3.75<=x<=3.75.

        real_t
        B3M_I0_PA_AS982(real_t x);    // |error| < 1.9e-07 for x>3.75.

        real_t
        B3M_I1_PA_AS983(real_t x);    // |error| < 8.0e-09 for -3.75<=x<=3.75.

        real_t
        B3M_I1_PA_AS984(real_t x);    // |error| < 2.2e-07 for x>3.75.

        real_t
        B_K0_PA_AS985(real_t x);      // |error| < 1.0e-08 for 0<x<=2.

        real_t
        B_K0_PA_AS986(real_t x);      // |error| < 1.9e-07 for x>2.

        real_t
        B_K1_PA_AS987(real_t x);      // |error| < 8.0e-09 for 0<x<=2.

        real_t
        B_K1_PA_AS988(real_t x);      // |error| < 2.2e-07 for x>2.

    }
}
#endif //UTIL_BESSEL_HPP
