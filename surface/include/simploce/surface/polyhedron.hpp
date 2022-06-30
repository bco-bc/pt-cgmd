/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#ifndef SURFACE_POLYHEDRON_HPP
#define SURFACE_POLYHEDRON_HPP

#include "./types/srf-types.hpp"
#include <cstdlib>
#include <iostream>

namespace simploce {

    /**
     * Represents a convex or simply connected (convex) three-dimensional shape (a surface)
     * in three dimensions with polygonal faces, edges and vertices.
     * @see <a href="https://en.wikipedia.org/wiki/Polyhedron">Wikipedia</a>
     */
    class Polyhedron {
    public:

        /**
         * Creates a polyhedron.
         * @param vertices Vertices.
         * @param faces Faces.
         * @param edges Edges.
         * @return Polyhedron.
         */
        static polyhedron_ptr_t create(const std::vector<vertex_ptr_t>& vertices,
                                       const std::vector<face_ptr_t>& faces,
                                       const std::vector<edge_ptr_t>& edges);

        /**
         * Constructor. The arguments should represent a simply connected shape. That is, its
         * Euler characteristics must be 2.
         * @param vertices Vertices.
         * @param faces Faces.
         * @param edges Edges without redundancy.
         */
        Polyhedron(std::vector<vertex_ptr_t> vertices,
                   std::vector<face_ptr_t> faces,
                   std::vector<edge_ptr_t> edges);

        // Not copyable.
        Polyhedron(const Polyhedron&) = delete;
        Polyhedron& operator = (const Polyhedron&) = delete;

        /**
         * Returns Euler characteristic Xi, defined as Xi=V-E+F, where V is the number of vertices,
         * E is the number of edges, and F is the number of faces.
         * @return Number. Returns 2 for any simply connected polyhedron.
         * @see <a href="https://en.wikipedia.org/wiki/Polyhedron#Characteristics">Wikipedia</a>
         */
        std::size_t eulerCharacteristic() const;

        /**
         * Perform a task with vertices, faces, and edges. The given task must exposes the
         * operator:
         * <code>
         * R operator () (const std::vector<vertex_ptr_t>& vertices,
         *                const std::vector<face_ptr_t>& faces,
         *                const std::vector<edge_ptr_t>& edges);
         * @tparam R Return type.
         * @tparam TASK Task type.
         * @param task Task to be performed. This may be a lambda expression.
         * @return Result.
         */
        template <typename R, typename TASK>
        R doWithAll(const TASK& task) const { return task(vertices_, faces_, edges_); }

        /**
         * Returns total area.
         * @return Area.
         */
        area_t area() const;

        /**
         * Writes edges to an output stream.
         * @param stream Stream.
         */
        void writeEdges(std::ostream& stream);

        /**
         * Returns number of faces.
         * @return Number of faces.
         */
        std::size_t numberOfFaces() const;

    private:

        friend class polyhedron_generator;
        friend class triangulation;

        void resetEdges_();

        void resetUnitVectorAtVertices_();

        std::vector<vertex_ptr_t> vertices_;
        std::vector<face_ptr_t> faces_;
        std::vector<edge_ptr_t> edges_;

    };

    /**
     * Writes polyhedron to output stream. Only vertices and faces are written into the stream.
     * @param stream Stream.
     * @param polyhedron Polyhedron.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const Polyhedron& polyhedron);
}

#endif //SURFACE_POLYHEDRON_HPP
