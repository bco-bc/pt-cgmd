/*
 * Author: André H. Juffer.
 * Created on 25/05/2022, 22:27.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_UTILS_HPP
#define BEM_UTILS_HPP

#include "simploce/types/u-types.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {
    namespace util {

        /**
         * Returns the inverse Debije length in a NaCl electrolyte solution
         * @tparam V Value type (float or double)
         * @param eps Relative dielectric permittivity constant. The default value
         * is at room temperature.
         * @param I Ionic strength, in mol/l. The default value is that of blood
         * plasma.
         * @param temperature Temperature in K. Default value corresponds to
         * the room temperature.
         * @return Inverse Debije length, in nm^-1.
         */
        template <typename V>
        V inverseDebijeLength(const molarity_t& I = 0.15,
                              const V eps = 78.2,
                              const temperature_t& temperature = 293.0) {
            auto IS = I() / units::mu<real_t>::l_to_nm3;  // mol/nm^3.
            auto E0 = units::mu<V>::E0;
            auto F = units::mu<V>::F;
            auto R = units::mu<V>::R;
            auto NA = units::si<V>::NA;
            return std::sqrt((2.0 * IS * F * F / NA) / (eps * E0 * R * temperature()));
        }
    }
}

#endif //BEM_UTILS_HPP
