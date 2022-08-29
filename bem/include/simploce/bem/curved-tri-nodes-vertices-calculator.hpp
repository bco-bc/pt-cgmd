/*
 * Author: André H. Juffer.
 * Created on 13/07/2022
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_CURVED_TRI_NODES_VERTICES_CALCULATOR_HPP
#define BEM_CURVED_TRI_NODES_VERTICES_CALCULATOR_HPP

#include "bem-calculator.hpp"

namespace simploce {

    /**
     * Assumes that the surface consist of curved triangles. Any function over any triangle is
     * numerically integrated and is consequently not assumed to be constant over a given triangle.
     * Vertices of triangles serve as nodes.
     */
    class CurvedTriNodesVerticesCalculator : public bem_calculator {
    public:

        /**
         * Returns an instance of this BEM calculator.
         * @param param Parameters.
         * @param surface Triangulated surface.
         * @return BEM calculator.
         */
        static bem_calc_ptr_t create(const param_ptr_t& param, const surface_ptr_t& surface);

        /**
         * Constructor.
         * @param param Parameters.
         * @param surface Triangulated surface.
         */
        CurvedTriNodesVerticesCalculator(param_ptr_t param, surface_ptr_t surface);

        void surfaceMatrix() override;

        void rightHandSide(std::vector<position_t>& positions, std::vector<charge_t>& charges) override;

        void solve() override;

        std::vector<el_pot> reactionPotentialSolute(const std::vector<position_t> &points) override;

        std::vector<el_pot> reactionPotentialSolvent(const std::vector<position_t> &points) override;

    private:

        param_ptr_t param_;
        surface_ptr_t surface_;

    };
}

#endif //BEM_CURVED_TRI_NODES_VERTICES_CALCULATOR_HPP
