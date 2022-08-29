/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#ifndef SURFACE_FACE_HPP
#define SURFACE_FACE_HPP

#include "types/srf-types.hpp"
#include <vector>
#include <map>

namespace simploce {

    /**
     * Surface element that is part of the boundary of an object. Holds vertices and edges.
     */
    class Face {
    public:

        virtual ~Face();

        // Not copyable.
        Face(const Face&) = delete;
        Face& operator = (const Face&) = delete;

        /**
         * Returns vertices indices.
         * @return Indices.
         */
        std::vector<int> indices() const;

        /**
         * Returns vertices.
         * @return Vertices.
         */
        std::vector<vertex_ptr_t> vertices() const;

        /**
         * Returns edges.
         * @return Edges.
         */
        std::vector<edge_ptr_t> edges() const;

        /**
         * Returns total area.
         * @return Area.
         */
        virtual area_t area() = 0;

        /**
         * Returns outward unit normal vector. Calculates normal from vertices' normals.
         * @return Unit normal vector.
         */
        virtual normal_t normal();

        /**
         * Returns center. Calculated as midpoint of face.
         * @return Center, associated outward unit normal vector.
         */
        virtual std::pair<position_t, normal_t> center();

    protected:

        Face(std::vector<vertex_ptr_t> vertices, std::vector<edge_ptr_t> edges);

        virtual void validate_() = 0;

    private:

        friend class Polyhedron;

        /**
         * Resets edges.
         * @param edges Edges.
         */
        virtual void resetEdges(const std::map<id_t, edge_ptr_t>& edges);

        std::vector<vertex_ptr_t> vertices_;
        std::vector<edge_ptr_t> edges_;
        std::vector<int> indices_;

        normal_t normal_{};
        std::pair<position_t, normal_t> center_{};
        bool hasNormal_{false};
        bool hasCenter_{false};

    };
}

#endif //SURFACE_FACE_HPP
