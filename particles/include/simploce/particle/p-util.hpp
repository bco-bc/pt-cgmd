/*
 * Author: André H. Juffer.
 * Created on 17/11/2021, 14:15.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef PARTICLES_P_UTIL_HPP
#define PARTICLES_P_UTIL_HPP

#include "p-types.hpp"

namespace simploce {
    namespace util {

        /**
         * Assigns velocity to the given particle compatible with the given temperature.
         * @param particle Particle.
         * @param temperature Temperature.
         */
        void assignVelocity(p_ptr_t &particle, const temperature_t& temperature);

    }
}

#endif //PARTICLES_P_UTIL_HPP
