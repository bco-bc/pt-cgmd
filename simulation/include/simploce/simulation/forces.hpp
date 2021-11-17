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
        Forces(box_ptr_t box, bc_ptr_t bc, ff_ptr_t forceField);

        energy_t nonBonded(const std::vector<p_ptr_t> &all,
                           const PairLists &pairLists);

        energy_t bonded(const std::vector<p_ptr_t> &all,
                        const std::vector<pg_ptr_t> &groups);

    private:

        box_ptr_t box_;
        bc_ptr_t bc_;
        ff_ptr_t forceField_;

    };

}

#endif //SIMULATION_FORCES_HPP
