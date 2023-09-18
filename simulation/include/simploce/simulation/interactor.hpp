/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 12:31 PM
 */

#ifndef INTERACTOR_HPP
#define INTERACTOR_HPP

#include "pair-list.hpp"
#include "simploce/types/s-types.hpp"

namespace simploce {

    /**
     * "One that interacts". Sits between the displacer and the forces calculator. It
     * ensures that the particle pair list is (re)evaluated on a regular basis.
     */
    class Interactor {
    public:

        /**
         * Constructor. All arguments are required.
         * @param param Simulation parameters.
         * @param forceField Force field.
         * @param PairListGenerator Particle pair list generator.
         * @param box Simulation box.
         */
        Interactor(param_ptr_t param,
                   pair_list_gen_ptr_t pairListGenerator,
                   forces_ptr_t forces);

        /**
         * Calculates forces on all particles in the given particle system.
         * @param simulationParameters Simulation parameters,
         * @param particleSystem Particle system.
         * @return Returns bonded, non-bonded, and external potential energy, respectively.
         */
        std::tuple<energy_t, energy_t, energy_t>
        interact(const p_system_ptr_t &particleSystem);

        /**
         * Returns the interaction energy of one particle with all other particles.
         * @param particle Particle.
         * @param particleSystem Particle system.
         * @return Returns bonded, non-bonded, and external interaction (potential) energy, respectively.
         */
        std::tuple<energy_t, energy_t, energy_t>
        interact(const p_ptr_t& particle,
                 const p_system_ptr_t &particleSystem);

        /**
         * Returns pair list.
         * @return Pair list.
         */
        pairlist_ptr_t
        pairList() const;

    private:

        param_ptr_t param_;
        pair_list_gen_ptr_t pairListGenerator_;
        forces_ptr_t forces_;

        pairlist_ptr_t pairList_;
    };


}

#endif /* INTERACTOR_HPP */

