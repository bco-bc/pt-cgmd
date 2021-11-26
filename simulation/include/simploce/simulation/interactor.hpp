/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 12:31 PM
 */

#ifndef INTERACTOR_HPP
#define INTERACTOR_HPP

#include "pair-lists.hpp"
#include "s-types.hpp"

namespace simploce {

    /**
     * "One that interacts". Its responsibility is ensure that the particle pair list
     * is (re)evaluated on a regular basis. Its sits in between a simulation and the
     * force calculation.
     */
    class Interactor {
    public:

        /**
         * Constructor. All arguments are required.
         * @param simulationParameters Simulation parameters.
         * @param forceField Force field.
         * @param PairListGenerator Particle pair list generator.
         * @param box Simulation box.
         */
        Interactor(sim_param_ptr_t simulationParameters,
                   pair_list_gen_ptr_t pairListGenerator,
                   forces_ptr_t forces);

        /**
         * Calculates forces on all particles in the given particle system.
         * @param simulationParameters Simulation parameters,
         * @param particleSystem Particle system.
         * @return Returns bonded and non-bonded potential energy, respectively.
         */
        std::pair<energy_t, energy_t> interact(const p_system_ptr_t &particleSystem);

        /**
         * Returns the interaction energy of one particle with all other particles.
         * @param particle Particle.
         * @param particleSystem Particle system.
         * @return Returns bonded and non-bonded interaction (potential) energy, respectively.
         */
        std::pair<energy_t, energy_t> interact(const p_ptr_t& particle,
                                               const p_system_ptr_t &particleSystem);

    private:

        sim_param_ptr_t simulationParameters_;
        pair_list_gen_ptr_t pairListsGenerator_;
        forces_ptr_t forces_;

        PairLists pairLists_;


    };


}

#endif /* INTERACTOR_HPP */

