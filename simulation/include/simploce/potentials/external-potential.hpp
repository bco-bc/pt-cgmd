/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#ifndef SIMULATION_EXTERNAL_POTENTIAL_HPP
#define SIMULATION_EXTERNAL_POTENTIAL_HPP

#include "../simulation/s-types.hpp"

namespace simploce {

    /**
     * An external force applied to -all- particles in the same way.
     */
    struct external_potential {

        virtual ~external_potential() = default;

        /**
         * Returns the potential energy and force on the particle due to this external
         * potential.
         * @param particle Particle.
         * @return Energy and force.
         */
        virtual std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) = 0;
    };
}

#endif //SIMULATION_EXTERNAL_POTENTIAL_HPP
