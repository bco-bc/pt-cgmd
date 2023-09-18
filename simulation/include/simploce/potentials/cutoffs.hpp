/*
 * Author: André H. Juffer.
 * Created on 24/08/2023, 13:46
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_CUTOFFS_HPP
#define SIMULATION_CUTOFFS_HPP

#include "simploce/types/s-types.hpp"

namespace simploce {

    /**
     * Holds cutoff distances for short range (SR) and long range (LR) interactions.
     */
    class Cutoffs {
    public:

        Cutoffs(dist_t sr, dist_t lr);

        // Not copyable.
        Cutoffs(const Cutoffs&) = delete;
        Cutoffs& operator = (const Cutoffs&) = delete;

        /**
         * Returns cutoff distance for short ranged interactions.
         * @return Cutoff distance.
         */
        dist_t shortRanged() const;

        /**
         * Returns cutoff distance for long ranged interactions.
         * @return Cutoff distance.
         */
        dist_t longRanged() const;

    private:

        dist_t sr_;
        dist_t lr_;
    };
}

#endif //SIMULATION_CUTOFFS_HPP
