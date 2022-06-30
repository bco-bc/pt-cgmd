/*
 * Author: André H. Juffer.
 * Created on 06/06/22, 21:30.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SURFACE_PENTAGON_HPP
#define SURFACE_PENTAGON_HPP

#include "face.hpp"
#include <vector>

namespace simploce {

    /**
     * A polygon with five edges and five vertices.
     */
    class Pentagon : public Face {
    public:

        /**
         * Create new pentagon. NOTE: New edges are constructed from the given vertices.
         * @param v1 Vertex.
         * @param v2 Vertex.
         * @param v3 Vertex.
         * @param v4 Vertex.
         * @param v5 Vertex.
         * @return Pentagon
         * @see <a href="https://en.wikipedia.org/wiki/Pentagon>Wikipedia</a>
         */
        static face_ptr_t create(const vertex_ptr_t& v1,
                                 const vertex_ptr_t& v2,
                                 const vertex_ptr_t& v3,
                                 const vertex_ptr_t& v4,
                                 const vertex_ptr_t& v5);

        /**
         * Constructor.
         * @param vertices Five vertices.
         * @param edges Five edges.
         */
        Pentagon(std::vector<vertex_ptr_t> vertices, std::vector<edge_ptr_t> edges);

        // Not copyable.
        Pentagon(const Pentagon&) = delete;
        Pentagon& operator = (const Pentagon&) = delete;

        /**
         * This will assume that the pentagon is flat, regular, and convex.
         * @return Area.
         * <a href="https://en.wikipedia.org/wiki/Pentagon#Derivation_of_the_area_formula">Wikipedia</a>
         */
        area_t area() const override;

        normal_t normal() const override;

    private:

        void validate_() override;

    };
}

#endif //SURFACE_PENTAGON_HPP
