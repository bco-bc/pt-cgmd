/*
 * Author: André H. Juffer.
 * Created on 15/07/2022
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_CURVE_HPP
#define BEM_CURVE_HPP

#include "types/bem-types.hpp"

namespace simploce {

    /**
     * Holds a parametrized curve between two vertices. A curve c(t) is given as c(t)=a+b*t+c*t^2+d*t^3, where
     * a, b, c and d are constant vectors, and t is a parameter ranging from 0 (starting point)
     * and 1 (end point) of the curve.
     */
    class Curve {
    public:

        /**
         * Creates a curve between two vertices.
         * @param start Start vertex.
         * @param end End vertex.
         * @return Curve.
         */
        static curve_ptr_t create(const vertex_ptr_t& start, const vertex_ptr_t& end);

        /**
         * Constructor
         * @param start Start vertex.
         * @param end End vertex.
         */
        Curve(vertex_ptr_t start, vertex_ptr_t end);

        /**
         * Returns a position with associated normal unit vector on this curve
         * @param t Parameter. Must be 0 <= t <= 1.
         * @return Position, unit vector.
         */
        std::pair<position_t, normal_t> point(real_t t);

    private:

        vertex_ptr_t start_;
        vertex_ptr_t end_;

        position_t a_{};
        position_t b_{};
        position_t c_{};
        position_t d_{};

    };
}

#endif //BEM_CURVE_HPP
