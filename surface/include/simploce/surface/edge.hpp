/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#ifndef SURFACE_EDGE_HPP
#define SURFACE_EDGE_HPP

#include "./types/srf-types.hpp"

namespace simploce {

    /**
     * Line segment joining two vertices.
     * @see https://en.wikipedia.org/wiki/Edge_(geometry)
     */
    class Edge {
    public:

        /**
         * Factory method.
         * @param start Start point.
         * @param end End point.
         */
        static edge_ptr_t create(const vertex_ptr_t& start, const vertex_ptr_t& end);

        /**
         * Constructor.
         * @param start Start point.
         * @param end End point.
         */
        Edge(vertex_ptr_t start, vertex_ptr_t end);

        // Noncopyable.
        Edge(const Edge&) = delete;
        Edge& operator = (const Edge&) = delete;

        /**
         * Returns an unique identifier for this edge.
         * @return Identifier.
         */
        id_t id() const;

        /**
         * Equality operator.
         * @param edge Edge.
         * @return Result.
         */
        bool operator == (const Edge& edge);

        /**
         * Returns start point.
         * @return Vertex.
         */
        vertex_ptr_t start() const;

        /**
         * Returns end point
         * @return Vertex.
         */
        vertex_ptr_t end() const;

        /**
         * Returns length of the edge.
         * @return Length.
         */
        length_t length() const;

    private:

        void validate_();

        id_t id_;
        vertex_ptr_t start_;
        vertex_ptr_t end_;
    };
}

#endif //SURFACE_EDGE_HPP
