/*
 * Author: André H. Juffer.
 * Created on 23/05/2022, 22:09.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/kernels.hpp"
#include "simploce/util/math-constants.hpp"

namespace simploce {
    namespace kernels {

        real_t Lij0(const real_t& epsRatio,
                    const position_t& r,
                    const normal_t& normal,
                    const position_t& r0) {
            const static auto FOUR_PI = 4.0 * simploce::math::constants<real_t>::PI;
            auto disv = r - r0;
            auto dis = norm<real_t>(disv);
            auto imp = inner<real_t>(disv, normal);
            auto f1 = 2.0 * (epsRatio - 1.0) / (epsRatio + 1.0);
            auto f2 = FOUR_PI * dis * dis * dis;
            return -f1 * imp / f2;
        }

        real_t Lij(const real_t& ka,
                   const real_t&epsRatio,
                   const position_t& r,
                   const normal_t &normal,
                   const position_t &r0) {
            const static auto PI = simploce::math::constants<real_t>::PI;
            auto disv = r - r0;
            auto dis = norm<real_t>(disv);
            auto imp = inner<real_t>(disv, normal);
            auto t0 = 1.0 + ka * dis;
            auto t1 = std::exp(-ka * dis);
            auto t2 = t0 * t1;
            auto t3 = 2.0 / (1.0 + epsRatio);
            auto t4 = 4.0 * PI * dis;
            auto t5 = t4 * dis * dis;
            return t3 * (1.0 - epsRatio * t2) * imp / t5;
        }

        std::tuple<real_t, real_t, real_t, real_t>
        KLMNij(const real_t& ka,
               const real_t& epsRatio,
               position_t &r,
               normal_t &n,
               position_t &r0,
               normal_t &n0) {
            const static auto PI = simploce::math::constants<real_t>::PI;
            auto disv = r - r0;                             // r-r0.
            auto dis = norm<real_t>(disv);                         // |r-r0|.
            auto dis2 = dis * dis;
            auto imp = inner<real_t>(disv, n);                 // (r-r0)*n.
            auto t0 = 1.0 + ka * dis;
            auto t1 = std::exp(-ka * dis);
            auto t2 = t0 * t1;
            auto t3 = 2.0 / (1.0 + epsRatio);
            auto t4 = 4.0 * PI * dis;
            auto t5 = t4 * dis2;
            auto t6 = t3 / t5;
            auto lij = t6 * (1.0 - epsRatio * t2) * imp;              // Kernel L(r,r0).
            auto kij = t3 * (1.0 - t1) / t4;                          // Kernel K(r,r0).
            auto imp0 = inner<real_t>(disv, n0);
            t3 *= epsRatio;
            auto nij = t6 * epsRatio * imp0 * (1.0 - t2 / epsRatio);  // Kernel N(r,r0).
            auto t7 = inner<real_t>(n, n0);
            auto t8 = imp * imp0;
            auto t9 = t8 / t5;
            auto m0 = 3.0 * t9 * (1.0 - t2) / dis2;
            auto m1 = t7 * (t2 - 1.0) / t5;
            auto m2 = t9 * ka * ka * t1;
            auto mij = t3 * (m0 + m1 -m2);                            // Kernel M(r;r0).
            return std::move(std::make_tuple(kij, lij, mij, nij));
        }
    }
}

