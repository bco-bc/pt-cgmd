/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_FORCES_HPP
#define SIMULATION_FORCES_HPP

#include "s-types.hpp"

namespace simploce {

    /**
     * Force calculator. Computes forces on particles.
     */
    class Forces {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @param forceField Force field.
         */
        Forces(bc_ptr_t bc, ff_ptr_t forceField);

        /**
         * Calculates forces on particles due to non-bonded interactions.
         * @param particleSystem Particle system.
         * @param pairLists Particle pair lists.
         * @return Potential energy.
         */
        energy_t nonBonded(const p_system_ptr_t& particleSystem,
                           const PairLists &pairLists);

        /**
         * Calculates forces on particles due to bonded interactions.
         * @param particleSystem Particle system.
         * @return Potential energy.
         */
        energy_t bonded(const p_system_ptr_t& particleSystem);

        /**
         * Returns interaction energy of one given particle with all other particles.
         * @param particle Particle.
         * @param particleSystem Particle system.
         * @return Bonded and non-bonded interaction (potential) energy.
         */
        std::pair<energy_t, energy_t> interaction(const p_ptr_t& particle,
                                                  const p_system_ptr_t& particleSystem);

    private:

        bc_ptr_t bc_;
        ff_ptr_t forceField_;

    };

}

#endif //SIMULATION_FORCES_HPP
