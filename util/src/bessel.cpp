/*
 * Author: André H. Juffer.
 * Created on 5/2/2024
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/bessel.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <cassert>
#include <cmath>
#include <limits>

namespace simploce {
    namespace math {

        static real_t epsilon = std::numeric_limits<real_t>::epsilon();

        real_t
        Bessel_J0(real_t x) {
            auto abs_x = std::fabs(x);
            if (abs_x < 8.0) {
                return boost::math::cyl_bessel_j(0, x);
                //return B1_Jn_S(0, x);             // Series expansion.
            } else {
                return boost::math::cyl_bessel_j(0, x);
                //return B1_J0_PA_AS943(x);               // |error| < 1.2e-08 for 3<=x<=8.
            }
        }

        real_t
        B1_Jn_S(size_t n, real_t &z) {
            assert(z < std::fabs(8.0) || "Must be -8.0 < z < 8.0");

            // Initialize.
            auto n_real = real_t(n);
            real_t nf = 1.0;        // nf = n! (Factorial).
            real_t two_pn = 1.0;     // 2^n
            real_t z_pow_n = 1.0;     // z^n
            for (size_t k = 1;  k <= n; ++k) {
                nf *= real_t(k);
                two_pn *= 2.0;
                z_pow_n *= z;
            }
            auto z2 = z * z;   // z * z
            auto z2_pow_k = 1.0;
            auto alpha = 1.0 / (nf * two_pn);    // For k=0;
            auto jn_k = alpha * z2_pow_k;       // For k=0;
            real_t k = 1.0;
            real_t jn_kP;
            do {
                jn_kP = jn_k;
                alpha *= -1.0/real_t(4.0 * k * (n_real + k));
                z2_pow_k *= z2;
                jn_k += (alpha * z2_pow_k);
                k += 1.0;
            } while (std::fabs(jn_k - jn_kP) > epsilon);
            jn_k *= z_pow_n;

            return jn_k;
        }

        real_t
        B1_J0_PA_AS941(real_t x) {
            assert (std::fabs(x) <= 3.0 || "Must be -3.0 <= x <= 3.0");

            auto a1 = x / 3.0;
            auto a2 = a1 * a1;
            auto a4 = a2 * a2;
            auto a6 = a4 * a2;
            auto a8 = a4 * a4;
            auto a10 = a8 * a2;
            auto a12 = a8 * a4;
            auto j0 = 1.0 - 2.2499997 * a2 + 1.2656208 * a4 - 0.3163866 * a6 + 0.0444479 * a8
                    - 0.0039444 * a10 + 0.0002100 * a12;

            return j0;
        }

        real_t
        B1_J0_PA_AS943(real_t x) {
            assert(x >= 3.0 || "Must be x >= 3.0.");

            auto a1 = 3.0 / x;
            auto a2 = a1 * a1;
            auto a3 = a2 * a1;
            auto a4 = a3 * a1;
            auto a5 = a4 * a1;
            auto a6 = a5 * a1;

            auto f0 = 0.79788456 - 0.00000077 * a1 - 0.00552740 * a2 - 0.00009512 * a3
                    + 0.00137237 * a4 - 0.00072805 * a5 + 0.00014476 * a6;
            auto theta = x - 0.78539816 - 0.04166397 * a1 - 0.00003954 * a2 + 0.00262573 * a3
                       - 0.00054125 * a4 - 0.00029333 * a5 + 0.00013558 * a6;
            auto j0 = f0 * std::cos(theta) / std::sqrt(x);

            return j0;
        }

        real_t
        Bessel_K0(real_t x) {
            assert(x >= epsilon || "Must be: x > epsilon");

            if (x > 700.0) return 0;
            if ( x <= 2.0)
                //return B_K0_PA_AS985(x);
                return boost::math::cyl_bessel_k(0, x);
            else
                //return B_K0_PA_AS986(x);
                return boost::math::cyl_bessel_k(0,x);
        }

        real_t
        Bessel_K1(real_t x) {
            assert(x >= epsilon || "Must be: x >= epsilon");

            if (x > 700.0) return 0;
            if ( x <= 2.0)
                //return B_K1_PA_AS987(x);
                return boost::math::cyl_bessel_k(1, x);
            else
                //return B_K1_PA_AS988(x);
                return boost::math::cyl_bessel_k(1, x);
        }

        real_t
        Bessel_I0(real_t x) {
            if (x <= 3.75)     // Between -3.75 and 3.75.
                return B3M_I0_PA_AS981(x);
            else               // Larger than 3.75.
                return B3M_I0_PA_AS982(x);
        }

        real_t
        Bessel_I1(real_t x) {
            if (x <= 3.75)     // Between -3.75 and 3.75.
                return B3M_I1_PA_AS983(x);
            else               // Larger than 3.75.
                return B3M_I1_PA_AS984(x);
        }

        real_t
        B3M_I0_PA_AS981(real_t x) {
            assert(std::fabs(x) <= 3.75 || "Must be -3.75 <= x <= 3.75");

            auto t = x / 3.75;
            auto t2 = t * t;
            auto t4 = t2 * t2;
            auto t6 = t4 * t2;
            auto t8 = t6 * t2;
            auto t10 = t8 * t2;
            auto t12 = t10 * t2;
            auto i0 = 1.0 + 3.5156229 * t2 + 3.0899424 * t4 + 1.2067492 * t6 + 0.2659732 * t8
                    + 0.0360768 * t10 + 0.0045813 * t12;

            return i0;
        }

        real_t
        B3M_I0_PA_AS982(real_t x) {
            assert (x > 3.75 || "Must be: x > 3.75.");

            auto b = std::sqrt(x);
            auto c = std::exp(-x);
            auto t = x / 3.75;
            auto a1 = 1 / t;
            auto a2 = a1 * a1;
            auto a3 = a2 * a1;
            auto a4 = a3 * a1;
            auto a5 = a4 * a1;
            auto a6 = a5 * a1;
            auto a7 = a6 * a1;
            auto a8 = a7 * a1;
            auto i0 = 0.39894228 + 0.01328592 * a1 + 0.00225319 * a2 - 0.00157565 * a3 + 0.00916281 * a4
                    - 0.02057706 * a5;
            i0 += 0.02635537 * a6 - 0.01647633 * a7 + 0.00392377 * a8;
            i0 *= 1.0 / (b * c);

            return i0;
        }

        real_t
        B3M_I1_PA_AS983(real_t x) {
            assert(fabs(x) <= 3.75 || "Must be: x <= 3.75");

            auto t =x / 3.75;
            auto t2 = t * t;
            auto t4 = t2 * t2;
            auto t6 = t4 * t2;
            auto t8 = t6 * t2;
            auto t10 = t8 * t2;
            auto t12 = t10 * t2;
            auto i1 = 0.5 + 0.87890594 * t2 + 0.51498869 * t4 + 0.15084934 * t6 + 0.02658733 * t8
                    + 0.00301532 * t10 + 0.00032411 * t12;
            i1 *= x;

            return i1;
        }

        real_t
        B3M_I1_PA_AS984(double x) {
            assert (x > 3.75 || "Must be x > 3.75");

            auto b = std::sqrt(x);
            auto c = std::exp(-x);
            auto t = x / 3.75;
            auto a1 =1/t;
            auto a2 = a1 * a1;
            auto a3 = a2 * a1;
            auto a4 = a3 * a1;
            auto a5 = a4 * a1;
            auto a6 = a5 * a1;
            auto a7 = a6 * a1;
            auto a8 = a7 * a1;
            auto i1 = 0.39894228 - 0.03988024 * a1 - 0.00362018 * a2 + 0.00163801 * a3 - 0.01031555 * a4;
            i1 += 0.02282967 * a5 - 0.02895312 * a6 + 0.01787654 * a7 - 0.00420059 * a8;
            i1 *= 1.0 /(b * c);

            return i1;
        }

        real_t
        B_K0_PA_AS985(real_t x) {
            assert(x >= epsilon && x <= 2.0 || "Must epsilon <= x <= 2.0");

            auto a = x / 2.0;
            auto b = std::log(x / 2.0);
            auto a2 = a * a;
            auto a4 = a2 * a2;
            auto a6 = a4 * a2;
            auto a8 = a6 * a2;
            auto a10 = a8 * a2;
            auto a12 = a10 * a2;
            auto i0 = Bessel_I0(x);
            auto k0 = -b * i0 - 0.57721566 + 0.42278420 * a2 + 0.23069756 * a4 + 0.03488590 * a6
                    + 0.00262698 * a8 + 0.00010750 * a10 + 0.00000740 * a12;

            return k0;
        }

        real_t
        B_K0_PA_AS986(real_t x) {
            assert(x > epsilon);

            auto b = std::sqrt(x);
            auto c = std::exp(x);
            auto a1 = 2.0 / x;
            auto a2 = a1 * a1;
            auto a3 = a2 * a1;
            auto a4 = a3 * a1;
            auto a5 = a4 * a1;
            auto a6 = a5 * a1;
            auto k0 = 1.25331414 - 0.07832358 * a1 + 0.02189568 * a2 - 0.01062446 * a3
                    + 0.00587872 * a4 - 0.00251540 * a5 + 0.00053208 * a6;
            k0 *= 1.0 / (c * b);

            return k0;
        }

        real_t
        B_K1_PA_AS987(real_t x) {
            assert(x >= epsilon && x <= 2.0 || "Must be: epsilon <= 2.0 <= 2.0");

            auto a = x / 2.0;
            auto b = std::log(x/2.0);
            auto a2 = a * a;
            auto a4 = a2 * a2;
            auto a6 = a4 * a2;
            auto a8 = a6 * a2;
            auto a10 = a8 * a2;
            auto a12 = a10 * a2;
            auto i1 = Bessel_I1(x);
            auto k1 = x * b * i1 + 1.0 + 0.15443144 * a2 - 0.67278579 * a4 - 0.18156897 * a6
                    - 0.01919402 * a8 - 0.00110404 * a10 - 0.00004686 * a12;
            k1 /= x;

            return k1;
        }

        real_t
        B_K1_PA_AS988(real_t x) {
            assert(x > epsilon || "Must be: x > epsilon");

            auto a1 = 2.0 / x;
            auto b = std::sqrt(x);
            auto c = std::exp(x);
            auto a2 = a1 * a1;
            auto a3 = a2 * a1;
            auto a4 = a3 * a1;
            auto a5 = a4 * a1;
            auto a6 = a5 * a1;
            auto k1 = 1.25331414 + 0.23498619 * a1 - 0.03655620 * a2 + 0.01504268 * a3
                    - 0.00780353 * a4 + 0.00325614 * a5 - 0.00068245 * a6;
            k1 *= 1.0 / (c * b);

            return k1;
        }

    }
}