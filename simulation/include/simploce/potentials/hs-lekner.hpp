/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 7 February 2024
 */

#ifndef SIMULATION_HS_LEKNER_HPP
#define SIMULATION_HS_LEKNER_HPP

#include "pair-potential.hpp"
#include "../conf/s-conf.hpp"

namespace simploce {

    /**
     * Pair interaction calculated according to the Lekner summation and a hard-sphere potential. Currently, only
     * available for Monte Carlo simulations. Force calculations were not implemented.
     */
    class HardSphereLekner : public pair_potential {
    public:

        HardSphereLekner(ff_ptr_t forceField,
                         bc_ptr_t bc,
                         lekner_ptr_t lekner);

        /*
         * @return Energy (forces are always zero).
         */
        std::pair<energy_t, force_t> operator () (const p_ptr_t &pi,
                                                  const p_ptr_t &pj) override;

        private:

            ff_ptr_t forceField_;
            bc_ptr_t bc_;
            lekner_ptr_t lekner_;
    };
}

#endif //SIMULATION_HS_LEKNER_HPP
