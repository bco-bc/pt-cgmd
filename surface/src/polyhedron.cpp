/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/surface/face.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/util/util.hpp"
#include "simploce/conf/u-conf.hpp"
#include <stdexcept>
#include <limits>
#include <map>

namespace simploce {

    polyhedron_ptr_t Polyhedron::create(const std::vector<vertex_ptr_t>& vertices,
                                        const std::vector<face_ptr_t>& faces,
                                        const std::vector<edge_ptr_t>& edges) {
        return std::make_shared<Polyhedron>(vertices, faces, edges);
    }

    Polyhedron::Polyhedron(std::vector<vertex_ptr_t> vertices,
                           std::vector<face_ptr_t> faces,
                           std::vector<edge_ptr_t> edges) :
        vertices_{std::move(vertices)}, faces_{std::move(faces)}, edges_{std::move(edges)} {
        util::Logger logger("simploce::Polyhedron::Polyhedron()");
        logger.trace("Entering.");

        auto xi = this->eulerCharacteristic();
        if (xi != 2) {
            util::logAndThrow(
            logger,
            "Polyhedron: Not a simply connected polyhedron: Euler characteristic must be 2.");
        }
        resetEdges_();

        logger.trace("Leaving.");
    }

    std::size_t
    Polyhedron::eulerCharacteristic() const {
        return vertices_.size() - edges_.size() + faces_.size();
    }

    area_t Polyhedron::area() const {
        area_t area{0.0};
        for (auto& face: faces_) {
            area += face->area();
        }
        return area;
    }

    void Polyhedron::writeEdges(std::ostream& stream) {
        stream << edges_.size() << std::endl;
        for (auto& edge: edges_) {
            stream << edge->start()->position() << edge->end()->position() << std::endl;
        }
    }

    std::size_t
    Polyhedron::numberOfFaces() const {
        return faces_.size();
    }

    void
    Polyhedron::resetEdges_() {
        util::Logger logger("simploce::Polyhedron::resetEdges_()");
        logger.trace("Entering");

        std::map<id_t, edge_ptr_t> reduced{};
        for (auto& edge: edges_) {
            auto pair = std::make_pair(edge->id(), edge);
            reduced.emplace(pair);
        }
        for (auto& face: faces_) {
            face->resetEdges(reduced);
        }

        logger.trace("Leaving.");
    }

    void
    Polyhedron::resetUnitVectorAtVertices_() {
        util::Logger logger("simploce::Polyhedron::resetUnitVectorAtVertices_()");
        logger.trace("Entering.");

        // Compute unit normal vector as an average of unit normal of faces.
        std::vector<normal_t> normals(vertices_.size(), normal_t{0.0, 0.0,0.0});
        std::vector<std::size_t> counters(vertices_.size(), 0);
        for (auto& face: faces_) {
            auto normal = face->normal();
            auto vertices = face->vertices();
            for (auto& vertex: vertices) {
                auto index = vertex->index();
                normals[index] += normal;
                counters[index] += 1;
            }
        }
        std::size_t counter = 0;
        for (auto& normal: normals) {
            normal /= counters[counter];
            counter += 1;
        }

        // Update vertices.
        for (auto& vertex: vertices_) {
            auto index = vertex->index();
            vertex->normal_(normals[index]);
        }

        logger.trace("Leaving.");
    }

    std::ostream& operator << (std::ostream& stream, const Polyhedron& polyhedron) {
        polyhedron.doWithAll<void>([&stream] (const std::vector<vertex_ptr_t>& vertices,
                                                   const std::vector<face_ptr_t>& faces,
                                                   const std::vector<edge_ptr_t>& edges) {
            stream << vertices.size() << conf::SPACE << faces.size() << std::endl;
            for (const auto& vertex: vertices) {
                stream << vertex->position() << std::endl;
            }
            for (const auto& vertex: vertices) {
                stream << vertex->normal() << std::endl;
            }
            for (const auto& face: faces) {
                for (const auto& vertex: face->vertices()) {
                    stream << conf::SPACE << vertex->index();
                }
                stream << std::endl;
            }
        });
        return stream;
    }
}