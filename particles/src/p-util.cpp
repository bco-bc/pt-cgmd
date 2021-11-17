/*
 * Author: André H. Juffer.
 * Created on 17/11/2021, 14:17.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/particle/p-util.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/util.hpp"
#include <random>

namespace simploce {
    namespace util {

        // From Maxwell velocity distribution.
        // https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
        void assignVelocity(p_ptr_t &particle, const temperature_t& temperature) {
            static real_t KB = units::mu<real_t>::KB;
            static std::random_device rd{};
            static std::mt19937 gen{rd()};
            static bool init = false;
            if ( !init ) {
                gen.seed(util::seedValue<std::size_t>());  // seed() From util.hpp
                init = true;
            }

            // For each velocity component, take its value from a Gaussian density for
            // velocity with standard deviation 'sigma' and zero average.
            real_t sigma = std::sqrt(KB * temperature / particle->mass());
            std::normal_distribution<real_t> gaussian{0.0, sigma};
            velocity_t v{};
            for (std::size_t k = 0; k != 3; ++k) {
                v[k] = gaussian(gen);
            }
            particle->velocity(v);
        }
    }
}

