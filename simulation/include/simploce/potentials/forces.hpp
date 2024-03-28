/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_FORCES_HPP
#define SIMULATION_FORCES_HPP

#include "simploce/types/s-types.hpp"

namespace simploce {

    /**
     * Force calculator. Computes forces on particles and associated interaction energies.
     */
    class Forces {
    public:

        /**
         * Constructor.
         * @param param Parameters
         * @param bc Boundary conditions.
         * @param forceField Force field.
         */
        Forces(param_ptr_t param,
               bc_ptr_t bc,
               ff_ptr_t forceField);

        /**
         * Calculates forces on particles due to non-bonded interactions.
         * @param particleSystem Particle system.
         * @param pairList Particle pair lists.
         * @return Potential energy.
         */
        energy_t
        nonBonded(const p_system_ptr_t& particleSystem,
                  const pairlist_ptr_t& pairList);

        /**
         * Calculates forces on particles to due external potentials.
         * @param particleSystem Particle system.
         * @return Potential energy.
         */
        energy_t
        external(const p_system_ptr_t& particleSystem);

        /**
         * Calculates forces on particles due to bonded interactions.
         * @param particleSystem Particle system.
         * @return Potential energy.
         */
        energy_t
        bonded(const p_system_ptr_t& particleSystem);

        /**
         * Returns the interaction energy of a single particle with all remaining particles in the particle
         * system.
         * @param particle Particle.
         * @param particleSystem Particle system.
         * @return Bonded, non-bonded, and external interaction (potential) energy.
         */
        std::tuple<energy_t, energy_t, energy_t>
        interaction(const p_ptr_t& particle,
                    const p_system_ptr_t& particleSystem);

        /**
         * Initiate. Always call this method first before updates.
         * @param particleSystem Particle system.
         */
        void
        initiate(const p_system_ptr_t& particleSystem);

        /**
         * Performs any update possibly required for calculating forces.
         * @param particleSystem Particle system
         */
        static void
        update(const p_system_ptr_t& particleSystem);

        /**
         * Performs any update possibly required for calculating forces.
         * @param particle Particle
         */
        static void
        update(const p_ptr_t& particle);

        /**
         * Fallback after updateStateAndAccumulated.
         */
        static void
        fallback();

        /**
         * Signals completion.
         */
        static void
        complete();

    private:

        param_ptr_t param_;
        rc_ptr_t cutoffs_;
        bc_ptr_t bc_;
        ff_ptr_t forceField_;
        bool mesoscopic_;
        bool concurrent_;

    };

}

#endif //SIMULATION_FORCES_HPP
