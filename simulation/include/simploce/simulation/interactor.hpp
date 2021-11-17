/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 12:31 PM
 */

#ifndef INTERACTOR_HPP
#define INTERACTOR_HPP

#include "pair-lists.hpp"
#include "forces.hpp"
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
                   ff_ptr_t forceField,
                   pair_list_gen_ptr_t pairListGenerator,
                   forces_ptr_t forces,
                   box_ptr_t box,
                   bc_ptr_t bc);

        /**
         * Calculates forces on all particles in the given particle system.
         * @param simulationParameters Simulation parameters,
         * @param particleSystem Particle system.
         * @return Returns non-bonded and bonded potential energy, respectively.
         */
        std::pair<energy_t, energy_t> interact(const p_system_ptr_t &particleSystem);

    private:

        sim_param_ptr_t simulationParameters_;
        ff_ptr_t forceField_;
        pair_list_gen_ptr_t pairListsGenerator_;
        forces_ptr_t forces_;
        box_ptr_t box_;
        bc_ptr_t bc_;

        PairLists pairLists_;


    };


}

#endif /* INTERACTOR_HPP */

