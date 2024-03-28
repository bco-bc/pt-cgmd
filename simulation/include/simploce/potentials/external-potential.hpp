/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#ifndef SIMULATION_EXTERNAL_POTENTIAL_HPP
#define SIMULATION_EXTERNAL_POTENTIAL_HPP

#include "simploce/types/s-types.hpp"

namespace simploce {

    /**
     * An external force, applied to -all- particles in the same way. External potentials may change due to
     * changes in the particle system and would require updates.
     */
    struct external_potential {

        virtual ~external_potential() = default;

        /**
         * Returns the potential energy and force on the particle due to this external potential.
         * @param particle Particle.
         * @return Energy and force. Note that the force is -not- assigned to particle, that is the
         * particle is -not- updated.
         */
        virtual std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) const = 0;

        /**
         * Initialize this potential, if required.
         * @param particleSystem Particle system.
         */
        virtual void
        initialize(const p_system_ptr_t& particleSystem) = 0;

        /**
         * Updates this potential, if required.
         * @param particleSystem Particle system.
         */
        virtual void
        update(const p_system_ptr_t& particleSystem) = 0;

        /**
         * Updates this potential, if required.
         * @param particle Particle.
         */
        virtual void
        update(const p_ptr_t& particle) = 0;

        /**
         * Falls back after an update.
         */
        virtual void
        fallback() = 0;

        /**
         * Completes this potential, e.g., after a simulation is completed.
         */
        virtual void
        complete() const = 0;
    };
}

#endif //SIMULATION_EXTERNAL_POTENTIAL_HPP
