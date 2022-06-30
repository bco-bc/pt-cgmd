/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/triangulation.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/surface/triangle.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/sphere.hpp"
#include "simploce/util/logger.hpp"
#include <vector>
#include <map>

namespace simploce {

    /**
     * Removes duplicated entries of edges in the given set of triangles.
     * @param triangles Triangles.
     * @return Collection of edges without redundancy.
     */
    static std::vector<edge_ptr_t>
    removeRedundancyInEdgesOfTriangles_(std::vector<face_ptr_t>& triangles) {
        using pair_t = std::pair<id_t, edge_ptr_t>;

        util::Logger logger("simploce::removeRedundancyInEdgesOfTriangles_()");
        logger.trace("Entering.");

        // Create unique set of edges.
        std::map<id_t, edge_ptr_t> current{};
        for (auto& triangle: triangles) {
            auto edges = triangle->edges();
            for (auto& edge: edges) {
                auto id = edge->id();
                pair_t pair = std::make_pair(id, edge);
                current.emplace(pair);
            }
        }

        // Place edges in a std::vector.
        std::vector<edge_ptr_t> reduced;
        for (const auto& pair : current) {
            reduced.emplace_back(pair.second);
        }
        logger.debug(std::to_string(reduced.size()) + ": Number of edges in non-redundant collection.");

        logger.trace("Leaving.");
        return std::move(reduced);
    }

    polyhedron_ptr_t
    triangulation::cubic(length_t sideLength) {
      real_t sl = sideLength();

      // 8 Vertices.
      std::vector<vertex_ptr_t> vertices;
      vertex_ptr_t v1{Vertex::create(position_t{0.0, 0.0, 0.0}, normal_t{})};
      vertex_ptr_t v2{Vertex::create(position_t{sl, 0.0, 0.0}, normal_t{})};
      vertex_ptr_t v3{Vertex::create(position_t{sl, sl, 0.0}, normal_t{})};
      vertex_ptr_t v4{Vertex::create(position_t {0.0, sl, 0.0}, normal_t{})};
      vertex_ptr_t v5{Vertex::create(position_t{0.0, 0.0, sl}, normal_t{})};
      vertex_ptr_t v6{Vertex::create(position_t {sl, 0.0, sl}, normal_t{})};
      vertex_ptr_t v7{Vertex::create(position_t {sl, sl, sl}, normal_t{})};
      vertex_ptr_t v8{Vertex::create(position_t {0.0, sl, sl}, normal_t{})};
      vertices.emplace_back(v1);
      vertices.emplace_back(v2);
      vertices.emplace_back(v3);
      vertices.emplace_back(v4);
      vertices.emplace_back(v5);
      vertices.emplace_back(v6);
      vertices.emplace_back(v7);
      vertices.emplace_back(v8);

       // 12 Triangles (each of the 6 faces of the cube in subdivided into two triangles).
      std::vector<face_ptr_t> triangles;
      auto t1 = Triangle::create(v1, v2, v3);
      auto t2 = Triangle::create(v1, v3, v4);
      auto t3 = Triangle::create(v3, v4, v8);
      auto t4 = Triangle::create(v3, v7, v8);
      auto t5 = Triangle::create(v1, v4, v8);
      auto  t6 = Triangle::create(v1, v5, v8);
      auto t7 = Triangle::create(v1, v2, v5);
      auto t8 = Triangle::create(v2, v5, v6);
      auto t9 = Triangle::create(v2, v3, v7);
      auto t10 = Triangle::create(v2, v6, v7);
      auto t11 = Triangle::create(v5, v6, v7);
      auto t12 = Triangle::create(v5, v7, v8);
      triangles.emplace_back(t1);
      triangles.emplace_back(t2);
      triangles.emplace_back(t3);
      triangles.emplace_back(t4);
      triangles.emplace_back(t5);
      triangles.emplace_back(t6);
      triangles.emplace_back(t7);
      triangles.emplace_back(t8);
      triangles.emplace_back(t9);
      triangles.emplace_back(t10);
      triangles.emplace_back(t11);
      triangles.emplace_back(t12);

      // Generate edges.
      auto edges = removeRedundancyInEdgesOfTriangles_(triangles);

      return std::make_shared<Polyhedron>(vertices, triangles, edges);
   }

    polyhedron_ptr_t
    triangulation::spherical(radius_t radius,  std::size_t numberOfTriangles) {
        util::Logger logger("simploce::triangulation::spherical()");
        logger.trace("Entering");

        // Create a dodecahedron.
        auto pair = sphere::dodecahedron(radius);
        auto& vertices = pair.first;
        auto& pentagons = pair.second;
        auto edges = removeRedundancyInEdgesOfTriangles_(pentagons);
        auto polyhedron = std::make_shared<Polyhedron>(vertices, pentagons, edges);
        auto area = polyhedron->area();
        logger.debug(std::to_string(area()) + ": Area dodecahedron.");
        logger.debug(std::to_string(polyhedron->eulerCharacteristic()) +
                     ": Euler characteristic of dodecahedron.");

        // Divide each pentagon into 5 triangles to generate 60 triangles and 32 vertices.
        pair = sphere::triangles60(polyhedron, radius);
        vertices = pair.first;
        auto triangles = pair.second;
        edges = removeRedundancyInEdgesOfTriangles_(triangles);
        polyhedron = std::make_shared<Polyhedron>(vertices, triangles, edges);
        area = polyhedron->area();
        logger.debug(std::to_string(area()) + ": Area triangulated sphere with 60 triangles.");

        // Continue dividing set of triangles until requested number of triangles is reached.
        pair = sphere::divide(polyhedron, numberOfTriangles, radius);
        vertices = pair.first;
        triangles = pair.second;

        // Edges.
        edges = removeRedundancyInEdgesOfTriangles_(triangles);
        auto surface = std::make_shared<Polyhedron>(vertices, triangles, edges);
        surface->resetUnitVectorAtVertices_();

        logger.trace("Leaving.");
        return std::move(surface);
    }

    polyhedron_ptr_t
    triangulation::parse(std::istream& stream) {
        int nVertices, nTriangles;
        std::string stringBuffer;

        stream >> nVertices >> nTriangles;
        std::getline(stream, stringBuffer);  // Read EOL.

        // Read vertices.
        // Positions.
        std::vector<position_t> positions;
        for (std::size_t k = 0; k != nVertices; ++k) {
            real_t x, y, z;
            stream >> x >> y >> z;
            std::getline(stream, stringBuffer);  // Read EOL.
            position_t r{x, y, z};
            positions.emplace_back(r);
        }
        // Normals.
        std::vector<normal_t> normals;
        for (std::size_t k = 0; k != nVertices; ++k) {
            real_t x, y, z;
            stream >> x >> y >> z;
            std::getline(stream, stringBuffer);  // Read EOL.
            normal_t un{x, y, z};
            normals.emplace_back(un);
        }
        // Create the vertices.
        std::vector<vertex_ptr_t> vertices;
        for (std::size_t k = 0; k != nVertices; ++k) {
            auto vertex = Vertex::create(positions[k], normals[k]);
            vertices.emplace_back(vertex);
        }

        // Read triangles.
        std::vector<face_ptr_t> triangles;
        for (std::size_t k = 0; k != nTriangles; ++k) {
            int i1, i2, i3;
            stream >> i1 >> i2 >> i3;
            std::getline(stream, stringBuffer);  // Read EOL.
            auto triangle = Triangle::create(vertices[i1], vertices[i2], vertices[i3]);
            triangles.emplace_back(triangle);
        }

        // Generate edges.
        auto edges = removeRedundancyInEdgesOfTriangles_(triangles);

        // Done
        return std::make_shared<Polyhedron>(vertices, triangles, edges);
    }
}
