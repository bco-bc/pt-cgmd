/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_LJ_HPP
#define SIMULATION_LJ_HPP

#include "pair-potential.hpp"
#include "s-types.hpp"
#include <utility>

namespace simploce {

    /**
     * Lennard Jones interaction.
     */
    class LJ : public pair_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param forceField Force field.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        LJ(ff_ptr_t forceField, box_ptr_t box, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        box_ptr_t box_;
        bc_ptr_t bc_;
    };


}

#endif //SIMULATION_LJ_HPP
