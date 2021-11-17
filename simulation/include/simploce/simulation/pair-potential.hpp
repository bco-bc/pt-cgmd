/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_PAIR_POTENTIAL_HPP
#define SIMULATION_PAIR_POTENTIAL_HPP

#include "s-types.hpp"

namespace simploce {

    /**
     * Interaction potential between two particles.
     * @tparam P Particle type.
     */
    struct pair_potential {

        /**
         * Returns potential energy and the force due to an interaction between two particles.
         * @param p1 Particle #1.
         * @param p2 Particle #2.
         * @return Potential energy and force on particle 1.
         */
        virtual std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) = 0;
    };
}

#endif //SIMULATION_PAIR_POTENTIAL_HPP
