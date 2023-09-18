/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 6/1/2022.
 */

#ifndef SIMULATION_HP_SR_HPP
#define SIMULATION_HP_SR_HPP

#include "pair-potential.hpp"
#include "soft-repulsion.hpp"
#include "hp.hpp"

namespace simploce {

    /**
     * Combined potential of the harmonic and soft repulsion pair potential.
     * This is a bonded potential for mesoscopic simulations, e.g., DPD.
     */
    class HarmonicSoftRepulsion : public pair_potential {
    public:

        /**
         * Constructor.
         * @param forceField Force field.
         * @param bc Boundary condition.
         * @param softRepulsion Soft repulsion potential.
         * @param harmonic Harmonic interaction potential.
         */
        HarmonicSoftRepulsion(ff_ptr_t forceField,
                              bc_ptr_t bc,
                              std::shared_ptr<SoftRepulsion> softRepulsion,
                              std::shared_ptr<HP> harmonic);

        std::pair<energy_t, force_t> operator () (const p_ptr_t &p1, const p_ptr_t &p2) override;

    private:

        ff_ptr_t forceField_;
        bc_ptr_t bc_;

        std::shared_ptr<SoftRepulsion> softRepulsion_;
        std::shared_ptr<HP> harmonic_;

    };
}

#endif //SIMULATION_HP_SR_HPP
