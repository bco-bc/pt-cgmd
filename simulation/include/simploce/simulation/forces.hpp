/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_FORCES_HPP
#define SIMULATION_FORCES_HPP

#include "s-types.hpp"
#include "pair-potential.hpp"

namespace simploce {

    /**
     * Computes forces on particles.
     */
    class Forces {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @param forceField Force field.
         */
        Forces(box_ptr_t box, bc_ptr_t bc, ff_ptr_t forceField);

        energy_t nonBonded(const std::vector<atom_ptr_t> &all,
                           const PairLists<Atom> &pairLists);

        energy_t bonded(const std::vector<atom_ptr_t> &all,
                        const std::vector<std::shared_ptr<ParticleGroup<Atom>>> &groups);

        energy_t nonBonded(const std::vector<bead_ptr_t> &all,
                           const PairLists<Bead> &pairLists);

        energy_t bonded(const std::vector<bead_ptr_t> &all,
                        const std::vector<std::shared_ptr<ParticleGroup<Bead>>> &groups);

    private:

        box_ptr_t box_;
        bc_ptr_t bc_;
        ff_ptr_t forceField_;

    };

}

#endif //SIMULATION_FORCES_HPP
