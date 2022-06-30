/*
 * Author: André H. Juffer.
 * Created on 06/06/22, 22:34.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/surface/sphere.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/types/srf-types.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/surface/pentagon.hpp"
#include "simploce/surface/triangle.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/map2.hpp"
#include <cmath>
#include <vector>

namespace simploce {
    namespace sphere {

        std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>>
        dodecahedron(const radius_t& radius) {
            util::Logger logger("simploce::sphere::dodecahedron()");
            logger.trace("Entering");

            // Compute side length (a) of each pentagon and of a cube (b).
            auto a= 2.0 * radius() / std::sqrt(3.0);
            auto b= 4.0 * radius() / ((1.0 + std::sqrt(5.0)) * std::sqrt(3.0));

            auto d = b * std::sqrt((3.0 - std::sqrt(5.0)) / 8.0);

            // 20 Vertices.
            std::vector<vertex_ptr_t> vertices;

            auto r = position_t{0.5 * a - d, 0.0, 0.5 * a + 0.5 * b};
            auto normal= normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v1{Vertex::create(r, normal)};
            vertices.emplace_back(v1);

            r = position_t{0.5 * a, 0.5 * a, 0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v2{Vertex::create(r, normal)};
            vertices.emplace_back(v2);

            r = position_t{0, 0.5 * a + 0.5 * b, 0.5 * a - d};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v3{Vertex::create(r,normal)};
            vertices.emplace_back(v3);

            r = position_t {-0.5 * a, 0.5 * a, 0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v4{Vertex::create(r, normal)};
            vertices.emplace_back(v4);

            r = position_t{-0.5 * a + d, 0.0, 0.5 * a + 0.5 * b};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v5{Vertex::create(r, normal)};
            vertices.emplace_back(v5);

            r = position_t {0.0, 0.5 * a + 0.5 * b, -0.5 * a + d};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v6{Vertex::create(r, normal)};
            vertices.emplace_back(v6);

            r = position_t{0.5 * a, 0.5 * a, -0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v7{Vertex::create(r, normal)};
            vertices.emplace_back(v7);

            r = position_t{0.5 * a - d, 0.0, -0.5 * a - 0.5 * b};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v8{Vertex::create(r, normal)};
            vertices.emplace_back(v8);

            r = position_t{-0.5 * a + d, 0.0, -0.5 * a - 0.5 * b};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v9{Vertex::create(r, normal)};
            vertices.emplace_back(v9);

            r = position_t{-0.5 * a, 0.5 * a, -0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v10{Vertex::create(r, normal)};
            vertices.emplace_back(v10);

            r = position_t{0.5 * a, -0.5 * a, -0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v11{Vertex::create(r, normal)};
            vertices.emplace_back(v11);

            r = position_t{-0.5 * a, -0.5 * a, -0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v12{Vertex::create(r, normal)};
            vertices.emplace_back(v12);

            r = position_t{0.0, -0.5 * a - 0.5 * b, -0.5 * a + d };
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v13{Vertex::create(r, normal)};
            vertices.emplace_back(v13);

            r = position_t{0.5 * a, -0.5 * a, 0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v14{Vertex::create(r, normal)};
            vertices.emplace_back(v14);

            r = position_t{0.0, -0.5 * a - 0.5 * b, 0.5 * a - d};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v15{Vertex::create(r, normal)};
            vertices.emplace_back(v15);

            r = position_t{-0.5 * a, -0.5 * a, 0.5 * a};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v16{Vertex::create(r, normal)};
            vertices.emplace_back(v16);

            r = position_t{0.5 * a + 0.5 * b, -0.5 * a + d, 0.0};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v17{Vertex::create(r, normal)};
            vertices.emplace_back(v17);

            r = position_t{0.5 * a + 0.5 * b, 0.5 * a - d, 0.0};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v18{Vertex::create(r, normal)};
            vertices.emplace_back(v18);

            r = position_t{-0.5 * a - 0.5 * b, -0.5 * a + d, 0.0};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v19{Vertex::create(r, normal)};
            vertices.emplace_back(v19);

            r = position_t{-0.5 * a - 0.5 * b, 0.5 * a - d, 0.0};
            normal = normal_t{r / norm<real_t>(r)};
            vertex_ptr_t v20{Vertex::create(r, normal)};
            vertices.emplace_back(v20);

            // 12 pentagons.
            std::vector<face_ptr_t> pentagons;
            auto p1 = Pentagon::create(v1, v2, v3, v4, v5);
            auto p2 = Pentagon::create(v5, v4, v20, v19, v16);
            auto p3 = Pentagon::create(v5, v16, v15, v14, v1);
            auto p4 = Pentagon::create(v1, v14, v17, v18, v2);
            auto p5 = Pentagon::create(v8, v7, v6, v10, v9);
            auto p6 = Pentagon::create(v9, v10, v20, v19, v12);
            auto p7 = Pentagon::create(v9, v12, v13, v11, v8);
            auto p8 = Pentagon::create(v8, v11, v17, v18, v7);
            auto p9 = Pentagon::create(v18, v7, v6, v3, v2);
            auto p10 = Pentagon::create(v6, v10, v20, v4, v3);
            auto p11 = Pentagon::create(v17, v14, v15, v13, v11);
            auto p12 = Pentagon::create(v13, v12, v19, v16, v15);
            pentagons.emplace_back(p1);
            pentagons.emplace_back(p2);
            pentagons.emplace_back(p3);
            pentagons.emplace_back(p4);
            pentagons.emplace_back(p5);
            pentagons.emplace_back(p6);
            pentagons.emplace_back(p7);
            pentagons.emplace_back(p8);
            pentagons.emplace_back(p9);
            pentagons.emplace_back(p10);
            pentagons.emplace_back(p11);
            pentagons.emplace_back(p12);

            logger.debug("Generated " + std::to_string(vertices.size()) + " vertices.");
            logger.debug("Generated " + std::to_string(pentagons.size()) + " pentagons.");

            logger.trace("Leaving");
            return std::move(std::make_pair(vertices, pentagons));
        }

        std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>>
        triangles60(const polyhedron_ptr_t& polyhedron, const radius_t& radius) {
            using result_t = std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>>;

            util::Logger logger("simploce::sphere::triangles60()");
            logger.trace("Entering.");

            auto result = polyhedron->doWithAll<result_t>([radius, &logger] (const std::vector<vertex_ptr_t>& oldVertices,
                                                            const std::vector<face_ptr_t>& pentagons,
                                                            const std::vector<edge_ptr_t>& edges) {
                logger.trace("Creating 60 triangles from 12 pentagons.");
                logger.debug(std::to_string(oldVertices.size()) + ": Initial number of vertices.");
                logger.debug(std::to_string(pentagons.size()) + ": Initial number of pentagons.");

                auto vertices = oldVertices;
                std::vector<face_ptr_t> triangles;
                for (const auto& face: pentagons) {
                    // Average of 5 oldVertices.
                    auto faceVertices = face->vertices();
                    position_t center{0.0, 0.0, 0.0};
                    for (auto& vertex: faceVertices) {
                        center += vertex->position();
                    }
                    center /= face->vertices().size();

                    // Put center on surface of sphere.
                    auto length = norm<real_t>(center);
                    center *= (radius() / length);

                    // Center becomes new vertex.
                    normal_t normal = center / radius();
                    auto vertex = Vertex::create(center, normal);
                    vertices.emplace_back(vertex);
                    logger.trace(std::to_string(vertices.size()) + ": Current number of vertices.");

                    // Create 5 new triangles.
                    auto t1 = Triangle::create(faceVertices[0], faceVertices[1], vertex);
                    auto t2 = Triangle::create(faceVertices[1], faceVertices[2], vertex);
                    auto t3 = Triangle::create(faceVertices[2], faceVertices[3], vertex);
                    auto t4 = Triangle::create(faceVertices[3], faceVertices[4], vertex);
                    auto t5 = Triangle::create(faceVertices[4], faceVertices[0], vertex);
                    triangles.emplace_back(t1);
                    triangles.emplace_back(t2);
                    triangles.emplace_back(t3);
                    triangles.emplace_back(t4);
                    triangles.emplace_back(t5);
                    logger.trace(std::to_string(triangles.size()) + ": Current number of triangles.");
                }
                logger.debug(std::to_string(vertices.size()) + ": Final number of oldVertices.");
                logger.debug(std::to_string(triangles.size()) + ": Final number of triangles.");
                logger.trace("Done: Created 60 triangles from 12 pentagons.");
                return std::make_pair(vertices, triangles);
            });

            logger.trace("Leaving.");
            return std::move(result);
        }

        std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>>
        divide(const polyhedron_ptr_t& polyhedron, std::size_t numberOfTriangles, const radius_t& radius) {
            util::Logger logger("simploce::sphere::divide()");
            logger.trace("Entering.");

            using result_t = std::pair<std::vector<vertex_ptr_t>, std::vector<face_ptr_t>>;
            using map2_t = MatrixMap<int, vertex_ptr_t>;

            // Results in 4 times more triangles.
            auto result = polyhedron->doWithAll<result_t>([numberOfTriangles, radius, &logger] (
                    const std::vector<vertex_ptr_t>& oldVertices,
                    const std::vector<face_ptr_t>& oldTriangles,
                    const std::vector<edge_ptr_t>& oldEdges) {

                logger.debug(std::to_string(oldVertices.size()) + ": Initial number of vertices.");
                logger.debug(std::to_string(oldTriangles.size()) + ": Initial number of triangles.");
                logger.debug(std::to_string(numberOfTriangles) + ": Requested number of triangles.");
                logger.debug(std::to_string(numberOfTriangles/2 + 2) + ": Requested number of vertices.");

                // map2 holds midpoint vertices of edges.
                map2_t map2{};
                std::vector<vertex_ptr_t> vertices = oldVertices;
                std::vector<face_ptr_t> triangles{};          // Will hold final set of triangles.
                auto current = oldTriangles;  // Set of triangles currently being divided.
                auto ntr = triangles.size();
                while (ntr < numberOfTriangles) {
                    triangles.clear();
                    for (const auto&triangle : current) {
                        auto edges = triangle->edges();
                        for (auto& edge: edges) {
                            // Calculate center of edge, if required.
                            auto start = edge->start();
                            auto startIndex = start->index();
                            auto end = edge->end();
                            auto endIndex = end->index();
                            if ( !map2.contains(startIndex, endIndex) ) {
                                // Not available yet. Create it.
                                auto mp = start->position() + 0.5 * (end->position() - start->position());

                                //logger.debug("|mp-start| = " + std::to_string(norm<real_t>(mp - start->position())));
                                //logger.debug("|mp-end| = " + std::to_string(norm<real_t>(mp - end->position())));

                                // Put it on the surface of the sphere.
                                auto length = norm<real_t>(mp);
                                mp *= (radius() / length);

                                // Create a new vertex.
                                normal_t normal = mp / radius();
                                auto vertex = Vertex::create(mp, normal);
                                vertices.emplace_back(vertex);

                                // Store new vertex for later reference.
                                map2.add(startIndex, endIndex, vertex);
                                map2.add(endIndex, startIndex, vertex);
                            }
                        }
                        logger.debug(std::to_string(vertices.size()) + ": Current number of vertices.");

                        // Create 4 new triangles.
                        auto indices = triangle->indices();  // Vertices indices.
                        auto c1 = map2.at(indices[0], indices[1]);
                        auto c2 = map2.at(indices[1], indices[2]);
                        auto c3 = map2.at(indices[2], indices[0]);
                        auto t1 = Triangle::create(c1, c2, c3);
                        auto t2 = Triangle::create(c1, c2, vertices[indices[1]]);
                        auto t3 = Triangle::create(c2, c3, vertices[indices[2]]);
                        auto t4 = Triangle::create(c3, c1, vertices[indices[0]]);
                        triangles.emplace_back(t1);
                        triangles.emplace_back(t2);
                        triangles.emplace_back(t3);
                        triangles.emplace_back(t4);
                        logger.debug(std::to_string(triangles.size()) + ": Current number of triangles.");
                    }
                    current = triangles;
                    ntr = triangles.size();
                }

                logger.debug(std::to_string(vertices.size()) + ": Final number of vertices.");
                logger.debug(std::to_string(triangles.size()) + ": Final number of triangles.");

                return std::make_pair(vertices, triangles);
            });


            logger.trace("Leaving.");
            return std::move(result);
        }

    }
}
