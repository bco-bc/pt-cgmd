/*
 * Author: André H. Juffer.
 * Created on 23/05/2022, 22:06.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_KERNELS_HPP
#define BEM_KERNELS_HPP

#include "./types/bem-types.hpp"
#include <tuple>

namespace simploce {
    namespace kernels {

        /**
         * Computes kernel L(r,r0), for zero ionic strength.
         * @param epsRatio Ratio of dielectric constants.
         * @param r Position on boundary.
         * @param normal Unit normal vector at r.
         * @param r0 Position on boundary.
         * @return Value.
         */
        real_t Lij0(const real_t& epsRatio,
                    const position_t& r,
                    const normal_t& normal,
                    const position_t& r0);

        /**
         * Computes kernel L(r,r0), for non-zero ionic strength.
         * @param ka Inverse Debye screening length.
         * @param epsRatio Ratio of dielectric constants.
         * @param r Position on boundary.
         * @param n Unit normal vector at r.
         * @param r0 Position on boundary, collocation point.
         * @return Value
         */
        real_t Lij(const real_t& ka,
                   const real_t&epsRatio,
                   const position_t& r,
                   const normal_t &normal,
                   const position_t &r0);

        /**
         * Returns all values for kernels K(r,r0), L(r,r0), M(r,r0), and N(r,r0) for .
         * @param ka Inverse Debye screening length.
         * @param epsRatio Ratio of dielectric constants.
         * @param r Position on boundary.
         * @param n Unit normal vector at r.
         * @param r0 Position on boundary, collocation point.
         * @param n0 Unit normal vector at r0.
         * @return K(r,r0), L(r,r0), M(r,r0), and N(r,r0).
         */
        std::tuple<real_t, real_t, real_t, real_t>
        KLMNij(const real_t& ka,
               const real_t& epsRatio,
               position_t &r,
               normal_t &n,
               position_t &r0,
               normal_t &n0);

        real_t Kij();
    }
}

#endif //BEM_KERNELS_HPP
