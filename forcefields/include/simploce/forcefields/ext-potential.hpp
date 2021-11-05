/*
 * File: potential.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 16, 2021., 1:00 PM.
 */

#ifndef EXT_POTENTIAL_HPP
#define EXT_POTENTIAL_HPP

#include "simploce/particle/p-types.hpp"
#include <utility>
#include <memory>

namespace simploce {

    /**
     * External potential.
     * @param P Particle type.
     * @see <a href="https://en.wikipedia.org/wiki/Potential">Wikipedia</a>
     */
    template <typename P>
    struct ext_potential {

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = std::shared_ptr<P>;

        /**
         * Amount of work (potential energy) required to move a particle from a reference point
         * to a given position.
         * @param p Particle
         * @return work (potential energy).
         */
        virtual energy_t energy(const particle_ptr_t& p) const = 0;

        /**
         * Sets the force experienced by a particle when placed at a given
         * position.
         * @param p Particle.
         * @return Work (potential energy).
         */
        virtual energy_t force(particle_ptr_t& p) const = 0;

    };
}

#endif //EXT_POTENTIAL_HPP
