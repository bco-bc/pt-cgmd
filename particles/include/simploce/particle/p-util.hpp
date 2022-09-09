/*
 * Author: André H. Juffer.
 * Created on 17/11/2021, 14:15.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef PARTICLES_P_UTIL_HPP
#define PARTICLES_P_UTIL_HPP

#include "p-types.hpp"
#include "simploce/util/util.hpp"

namespace simploce {
    namespace util {

        /**
         * Assigns velocity to the given particle compatible with the given temperature.
         * @param particle Particle.
         * @param temperature Temperature.
         * @param isMesoscale If true, temperature is expressed as kT=n, where k is the
         * Boltzmann constant and n is a nonnegative integer.
         */
        void assignVelocity(p_ptr_t &particle, const temperature_t& temperature, bool isMesoscale = false);

        /**
         * Finds particle in a given particle collection. This is slow for large collections of particles when
         * repeatable used.
         * @param id Particle identifier.
         * @param particles Particles.
         * @return Particle or nullptr if the particle cannot be identified.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        std::shared_ptr<P> find(simploce::id_t id, const CONT<std::shared_ptr<P>>& particles)
        {
            for (auto& p : particles) {
                if ( p->id() == id ) {
                    return p;
                }
            }
            return nullptr;
        }

        /**
         * Ensure that the total linear moment of the given particle system adds up to zero.
         * @param particleSystem Particle system.
         */
        void removeOverallLinearMomentum(const p_system_ptr_t& particleSystem);
    }
}

#endif //PARTICLES_P_UTIL_HPP
