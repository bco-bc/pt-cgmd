/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 20:58.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/flat-tri-nodes-vertices-calculator.hpp"
#include "simploce/bem/bem-data.hpp"
#include "simploce/bem/kernels.hpp"
#include "simploce/bem/rhs.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/surface/face.hpp"
#include "simploce/util/logger.hpp"
#include <vector>
#include <stdexcept>
#include <memory>

namespace simploce {
    namespace bem_flat_tri_nodes_vertices {

        static std::shared_ptr<BEMData> bemData_{};

        /**
         * Assign vertices as nodes.
         * @param surface Surface.
         * @return Nodes.
         */
        static std::vector<BEMData::node_t>
        assignNodes_(const surface_ptr_t& surface) {
            auto nodes = surface->doWithAll<std::vector<BEMData::node_t>>([](
                    const std::vector<vertex_ptr_t> &vertices,
                    const std::vector<face_ptr_t> &faces,
                    const std::vector<edge_ptr_t> &edges) {
                std::vector<BEMData::node_t> nodes;
                for (auto& vertex: vertices) {
                    BEMData::node_t node;
                    node.r = vertex->position();
                    node.nv = vertex->normal();
                    node.index = vertex->index();
                    nodes.emplace_back(node);
                }
                return std::move(nodes);
            });
            return std::move(nodes);
        }

        static void
        initialize_(const param_ptr_t& param, const surface_ptr_t& surface) {
            util::Logger logger("simploce::flat_tri_nodes_vertices::initialize_()");
            logger.trace("Entering.");

            auto numberOfVertices = surface->numberOfVertices();
            bemData_ = std::make_shared<BEMData>(param, numberOfVertices);
            auto& data = *bemData_;
            data.nodes = assignNodes_(surface);

            logger.debug(std::to_string(surface->numberOfVertices()) + ": Number of vertices.");
            logger.debug(std::to_string(data.S.size()) + ": Number of elements surface matrix.");

            logger.trace("Leaving.");
        }

        static std::vector<el_pot>
        reactionPotential_(real_t factor,
                           const std::vector<position_t> &points,
                           const surface_ptr_t& surface) {
            util::Logger logger("simploce::bem_flat_tri_nodes_vertices::reactionPotential()");
            logger.trace("Entering.");

            logger.debug(std::to_string(points.size()) +
                         ": Number of points for electric potential calculation.");

            auto result =
                surface->doWithAll<std::vector<el_pot>>([points, factor](
                        const std::vector<vertex_ptr_t> &vertices,
                        const std::vector<face_ptr_t> &faces,
                        const std::vector<edge_ptr_t> &edges) {
                    auto& data = *bem_flat_tri_nodes_vertices::bemData_;
                    auto epsRatio = data.epsRatio;
                    std::vector<el_pot> rp(points.size(), 0.0);
                    for (auto k = 0; k != points.size(); ++k) {
                        auto rk = points[k];
                        for (const auto &face: faces) {
                            auto indices = face->indices();
                            auto area = face->area();
                            auto pair = face->center();
                            auto value = kernels::Lij0(epsRatio, pair.first, pair.second, rk);
                            value *= area() / 3.0;
                            for (int index: indices) {
                                auto x = data.x[index];
                                rp[k] += factor * value * x;
                            }
                        }
                    }
                    return std::move(rp);
                });

            logger.trace("Leaving.");
            return std::move(result);
        }

    }

    bem_calc_ptr_t
    FlatTriNodesVerticesCalculator::create(const param_ptr_t &param, const surface_ptr_t &surface) {
        return std::make_shared<FlatTriNodesVerticesCalculator>(param, surface);
    }

    FlatTriNodesVerticesCalculator::FlatTriNodesVerticesCalculator(param_ptr_t param, surface_ptr_t surface) :
            param_(std::move(param)), surface_{std::move(surface)} {
        util::Logger logger("simploce::FlatTriNodesVerticesCalculator::FlatTriNodesVerticesCalculator()");
        logger.trace("Entering");

        if (!param_) {
            util::logAndThrow(logger, "Missing parameters.");
        }
        if (!surface_) {
            util::logAndThrow(logger, "Missing (triangulated) surface.");
        }
        bem_flat_tri_nodes_vertices::initialize_(param_, surface_);

        logger.trace("Leaving.");
    }

    void
    FlatTriNodesVerticesCalculator::surfaceMatrix() {
        util::Logger logger("simploce::FlatTriNodesVerticesCalculator::surfaceMatrix()");
        logger.trace("Entering");

        surface_->doWithAll<void>([&logger] (const std::vector<vertex_ptr_t>& vertices,
                                                  const std::vector<face_ptr_t>& faces,
                                                  const std::vector<edge_ptr_t>& edges) {

            auto& data = *bem_flat_tri_nodes_vertices::bemData_;
            auto epsRatio = data.epsRatio;
            auto& S = data.S;
            auto size = S.size();
            logger.debug(std::to_string(size) + ": Size (number of elements) of surface matrix.");

            // Compute elements of S.
            S.setZero();
            for (auto i = 0; i != vertices.size(); ++i) {
                S(i,i) = 1.0;
            }
            auto& nodes = data.nodes;

            // No singularities to worry about.
            for (const auto& face: faces) {
                auto indices = face->indices();
                auto area = face->area();
                auto pair = face->center();  // Position center and outward normal vector.
                for (const auto& node_i: nodes) {
                    auto& index_i = node_i.index;
                    auto& r0 = node_i.r;
                    auto value = kernels::Lij0(epsRatio, pair.first, pair.second, r0);
                    value *= area() / 3.0;
                    for (auto index_j : indices) {
                        S(index_i, index_j) -= value;
                    }
                }
            }
            //std::clog << "Before LU: " << std::endl << S << std::endl;

            // LU-decomposition.
            data.lu.compute(S);

            //std::clog << "After LU: " << std::endl << S << std::endl;
        });

        // Done.
        logger.trace("Leaving.");
    }

    void
    FlatTriNodesVerticesCalculator::rightHandSide(std::vector<position_t>& positions, std::vector<charge_t>& charges) {
        static util::Logger logger("simploce::FlatTriNodesVerticesCalculator::rightHandSide");
        logger.trace("Entering");

        if (positions.size() != charges.size()) {
            util::logAndThrow(logger, "Length charges and positions vector are not the same.");
        }
        if (positions.empty() || charges.empty()) {
            util::logAndThrow(logger, "Positions and/or charges not provided.");
        }

        auto data = bem_flat_tri_nodes_vertices::bemData_;
        data->b = rhs::rightHandSide(data->nodes,
                                     positions,
                                     charges,
                                     data->epsSolvent,
                                     data->epsSolute);

        logger.trace("Leaving");
    }

    void
    FlatTriNodesVerticesCalculator::solve() {
        static util::Logger logger("simploce::FlatTriNodesVerticesCalculator::solve()");
        logger.trace("Entering.");
        auto& data = *bem_flat_tri_nodes_vertices::bemData_;
        data.x = data.lu.solve(data.b);
        //std::clog << "Solution: "<< std::endl << data.x << std::endl;
        logger.trace("Leaving.");
    }

    std::vector<el_pot>
    FlatTriNodesVerticesCalculator::reactionPotentialSolute(const std::vector<position_t> &points) {
        util::Logger logger("simploce::FlatTriNodesVerticesCalculator::reactionPotentialSolute()");
        logger.trace("Entering.");

        auto& data = *bem_flat_tri_nodes_vertices::bemData_;
        auto epsRatio = data.epsRatio;
        auto factor = (epsRatio + 1) / 2.0;
        auto result = bem_flat_tri_nodes_vertices::reactionPotential_(factor, points, surface_);

        logger.trace("Leaving.");
        return std::move(result);
    }

    std::vector<el_pot>
    FlatTriNodesVerticesCalculator::reactionPotentialSolvent(const std::vector<position_t> &points) {
        util::Logger logger("simploce::FlatTriNodesVerticesCalculator::reactionPotentialSolvent()");
        logger.trace("Entering.");

        auto& data = *bem_flat_tri_nodes_vertices::bemData_;
        auto epsRatio = data.epsRatio;
        auto factor = (epsRatio + 1) / (2.0 * epsRatio);
        auto result = bem_flat_tri_nodes_vertices::reactionPotential_(factor, points, surface_);

        logger.trace("Leaving.");
        return std::move(result);

    }

}