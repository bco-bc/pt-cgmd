/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/11/21.
 */

#ifndef SIMULATION_LJ_ELECTROSTATICS_HPP
#define SIMULATION_LJ_ELECTROSTATICS_HPP

#include "s-types.hpp"

namespace simploce {
    namespace forces {

        /**
         * Computes the Lennard-Jones and electrostatic interaction between the
         * atom pairs
         * @param pairLists Atom pair lists
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Potential energy.
         */
        energy_t computeLjElectrostatic(const std::vector<atom_ptr_t> &all,
                                        const PairLists<Atom> &pairLists,
                                        const ff_ptr_t &forceField,
                                        const box_ptr_t &box,
                                        const bc_ptr_t &bc);

        /**
         * Computes the Lennard-Jones and electrostatic interaction between the
         * bead pairs
         * @param pairLists Bead pair lists
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Potential energy.
         */
        energy_t computeLjElectrostatic(const std::vector<bead_ptr_t> &all,
                                        const PairLists<Bead> &pairLists,
                                        const ff_ptr_t &forceField,
                                        const box_ptr_t &box,
                                        const bc_ptr_t &bc);
    }

}
#endif //SIMULATION_LJ_ELECTROSTATICS_HPP
