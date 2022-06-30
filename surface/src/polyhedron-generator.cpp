/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/13/22.
 */

#include "simploce/surface/polyhedron-generator.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    polyhedron_generator::~polyhedron_generator() = default;

    void
    polyhedron_generator::mapOnto(std::vector<dot_t>& dots,
                                  polyhedron_ptr_t& polyhedron) {
        util::Logger logger("polyhedron_generator::mapOnto(const std::vector<dot_t>& dots)");
        logger.trace("Entering.");

        polyhedron->doWithAll<void>([&dots, &logger] (const std::vector<vertex_ptr_t>& vertices,
                                                           const std::vector<face_ptr_t>& faces,
                                                           const std::vector<edge_ptr_t>& edges) {

            // Move geometric center dots to origin.
            position_t dotsCenter{0.0, 0.0, 0.0};
            for (auto& dot: dots) {
                dotsCenter += dot;
            }
            dotsCenter /= dots.size();
            logger.debug(util::toString(dotsCenter) + ": Current geometric center of dots.");
            real_t largest{0.0};
            for (auto& dot: dots) {
                dot -= dotsCenter;
                auto length = norm<real_t>(dot);
                if (length > largest) {
                    largest = length;
                }
            }
            dotsCenter = {0.0, 0.0, 0.0};
            for (auto& dot: dots) {
                dotsCenter += dot;
            }
            dotsCenter /= dots.size();
            logger.debug(util::toString(dotsCenter) + ": New geometric center of dots.");
            logger.debug(std::to_string(largest) + ": Length longest dot.");

            // Move geometric center vertices to origin.
            position_t verticesCenter{0.0, 0.0, 0.0};
            for (const auto& v: vertices) {
                verticesCenter += v->position();
            }
            verticesCenter /= vertices.size();
            logger.debug(util::toString(verticesCenter) + ": Current geometric center of vertices.");
            for (auto& vertex: vertices) {
                auto r = vertex->position();
                r -=verticesCenter;
                vertex->position_(r);
            }
            verticesCenter = {0.0, 0.0, 0.0};
            for (const auto& v: vertices) {
                verticesCenter += v->position();
            }
            verticesCenter /= vertices.size();
            logger.debug(util::toString(verticesCenter) + ": New geometric center of vertices.");

            // Adjust length of vertices to dimension of dotted surface.
            for (const auto& vertex: vertices) {
                auto r = vertex->position();
                r *= (largest / norm<real_t>(r));
                vertex->position_(r);
            }

            // Find the 'best' dot to replace vertex.
            auto large2 = std::numeric_limits<real_t>::max();
            auto points = dots;
            for (const auto& vertex: vertices) {

                bool found{false};
                auto selected = points.begin();

                // Current vertex.
                auto r = vertex->position();
                auto vertexLength = norm<real_t>(r);
                bool sameQuadrant{false};                       // Are there dots in the same quadrant as given vertex?

                real_t closest2 = large2;
                for (auto iter = points.begin(); iter != points.end(); ++iter) {
                    auto dot = *iter;
                    if ( util::sgn(dot[0]) == util::sgn(r[0]) &&
                         util::sgn(dot[1]) == util::sgn(r[1]) &&
                         util::sgn(dot[2]) == util::sgn(r[2]) ) {
                        // Dot and vertex are in the same quadrant.
                        sameQuadrant = true;
                        auto dotLength = norm<real_t>(dot);
                        auto rv = (dotLength / vertexLength) * r;
                        auto distance2 = norm_square<real_t>(rv - dot);
                        if (distance2 < closest2) {
                            found = true;
                            closest2 = distance2;
                            selected = iter;
                        }
                    }
                }
                if (!found) {
                    if ( !sameQuadrant ) {
                        logger.debug("No dots in same quadrant as vertex.");
                    }
                    util::logAndThrow(logger,
                                         "Polyhedron: Could not map vertex onto dotted surface.");
                }
                auto point = *selected;
                vertex->position_(point);

                // Remove selected dot to avoid being selected again.
                points.erase(selected);
            }
        });
        polyhedron->resetUnitVectorAtVertices_();

        logger.trace("Leaving");
    }

}