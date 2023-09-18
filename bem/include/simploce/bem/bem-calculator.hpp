/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 14:27.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_BEM_CALCULATOR_HPP
#define BEM_BEM_CALCULATOR_HPP

#include "./types/bem-types.hpp"
#include <vector>

namespace simploce {

    /**
     * BEM Strategy for calculating electric potentials.
     */
    struct bem_calculator {

        virtual ~bem_calculator() {}

        /**
         * Computes the surface matrix S in Sx=b, where x represents the
         * set of unknowns at the collocation points and b is the source vector.
         */
        virtual void surfaceMatrix() = 0;

        /**
         * Computes right-hand-side b in Sx=b.
         * @param positions Positions of charges.
         * @param charges Charge values.
         */
        virtual void rightHandSide(std::vector<position_t>& positions, std::vector<charge_t>& charges) = 0;

        /**
         * Finds the solution x of Sx=b.
         */
        virtual void solve() = 0;

        /**
         * Returns reaction potential at specified points located in the solute region.
         * @param points Points.
         * @return Reaction potentials.
         */
        virtual std::vector<el_pot_t> reactionPotentialSolute(const std::vector<position_t> &points) = 0;

         /**
         * Returns reaction potential at specified points located in the solvent region.
         * @param points Points.
         * @return Reaction potentials.
         */
       virtual std::vector<el_pot_t> reactionPotentialSolvent(const std::vector<position_t> &points) = 0;
    };
}

#endif //BEM_BEM_CALCULATOR_HPP
