/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/12/21.
 */

#ifndef SIMULATION_LJ_HPP
#define SIMULATION_LJ_HPP

#include "simploce/potentials/pair-potential.hpp"
#include "simploce/simulation/s-types.hpp"
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
        LJ(ff_ptr_t forceField, bc_ptr_t bc);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        friend class LJ_RF;
        friend class LJ_SF;

        static std::pair<energy_t, force_t> forceAndEnergy(const dist_vect_t& rij,
                                                           real_t Rij,
                                                           real_t Rij2,
                                                           real_t C12,
                                                           real_t C6);

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
    };


}

#endif //SIMULATION_LJ_HPP
