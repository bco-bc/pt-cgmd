/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */


#ifndef SURFACE_TRIANGLE_HPP
#define SURFACE_TRIANGLE_HPP

#include "face.hpp"

namespace simploce {

    /**
     * A polygon with three edges and three vertices.
     */
    class Triangle: public Face {
    public:

        /**
         * Creates new triangle. NOTE: New edges are constructed from the given vertices.
         * @param v1 Vertex.
         * @param v2 Vertex.
         * @param v3 Vertex.
         * @return Triangle.
         */
        static face_ptr_t create(const vertex_ptr_t& v1,
                                 const vertex_ptr_t& v2,
                                 const vertex_ptr_t& v3);

        /**
         * Constructor.
         * @param vertices Three vertices.
         * @param edges Three edges.
         */
        Triangle(std::vector<vertex_ptr_t> vertices, std::vector<edge_ptr_t> edges);

        // Not copyable.
        Triangle(const Triangle&) = delete;
        Triangle& operator = (const Triangle&) = delete;

        area_t area() const override;

        normal_t normal() const override;

    private:

        void validate_() override;

    };
}

#endif //SURFACE_TRIANGLE_HPP
