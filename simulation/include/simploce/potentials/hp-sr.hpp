/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 6/1/2022.
 */

#ifndef SIMULATION_HP_SR_HPP
#define SIMULATION_HP_SR_HPP

#include "pair-potential.hpp"

namespace simploce {

    /**
     * Combination of the harmonic and soft repulsion pair potential. This is a bonded potential.
     */
    class HarmonicSoftRepulsion : public pair_potential {
    public:

        /**
         * Constructor.
         * @param forceField Force field.
         * @param bc Boundary condition.
         * @param cutoff Cutoff distance.
         */
        HarmonicSoftRepulsion(ff_ptr_t forceField, bc_ptr_t bc, dist_t cutoff);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;
        dist_t cutoff_;

    };
}

#endif //SIMULATION_HP_SR_HPP
