/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_PAIR_POTENTIAL_HPP
#define SIMULATION_PAIR_POTENTIAL_HPP

#include "simploce/simulation/s-types.hpp"

namespace simploce {

    /**
     * Interaction potential between two particles.
     * @tparam P Particle typeName.
     */
    struct pair_potential {

        virtual ~pair_potential() = default;

        /**
         * Returns potential energy and the force due to an interaction between two particles.
         * @param pi Particle.
         * @param pj Particle.
         * @return Potential energy and force on particle i. -Not- assigned.
         */
        virtual std::pair<energy_t, force_t> operator () (const p_ptr_t &pi, const p_ptr_t &pj) = 0;
    };
}

#endif //SIMULATION_PAIR_POTENTIAL_HPP
