/*
 * File: pair_potential.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 16, 2021., 4:56 PM.
 */

#ifndef PAIR_POTENTIAL_HPP
#define PAIR_POTENTIAL_HPP

#include "simploce/types/u-types.hpp"
#include <utility>
#include <memory>

namespace simploce {

    /**
     * Potential for two interacting particles.
     */
    template<typename P>
    struct pair_potential {

        /**
         * Particle pointer type.
         */
        using particle_ptr_t = std::shared_ptr<P>;

        /**
         * Returns amount of work required to place two interacting particles at a specified
         * configuration.
         * @return Interaction energy.
         */
         virtual energy_t energy(const particle_ptr_t& p1, const particle_ptr_t & p2) const = 0;

        /**
         * Sets the forces on two interacting particles.
         * @return Interaction energy in kJ/mol.
         * @param p1 Particle #1.
         * @param p2 Particle #2.
         * @return Interaction energy.
         */
        virtual energy_t force(particle_ptr_t& p1, particle_ptr_t& p2) const = 0;

    };

}

#endif //PAIR_POTENTIAL_HPP
