/*
 * Author: André H. Juffer.
 * Created on 13/07/2022
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/curved-tri-nodes-vertices-calculator.hpp"
#include "simploce/bem/bem-data.hpp"
#include "simploce/bem/rhs.hpp"
#include "simploce/bem/curve.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/surface/face.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/util/logger.hpp"
#include <memory>
#include <map>

namespace simploce {
    namespace bem_curved_tri_nodes_vertices {

        static std::shared_ptr<BEMData> bemData_{};

        // ID is the edge identifier.
        static std::map<id_t, curve_ptr_t> curves_{};

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
            util::Logger logger("simploce::bem_curved_tri_nodes_vertices::initialize_()");
            logger.trace("Entering.");

            logger.debug(std::to_string(surface->numberOfVertices()) + ": Number of vertices.");

            auto numberOfVertices = surface->numberOfVertices();
            bemData_ = std::make_shared<BEMData>(param, numberOfVertices);
            auto& data = *bemData_;

            data.nodes = assignNodes_(surface);
            surface->doWithAll<void>([&logger] (const std::vector<vertex_ptr_t>& vertices,
                                                     const std::vector<face_ptr_t>& faces,
                                                     const std::vector<edge_ptr_t>& edges) {
                logger.debug(std::to_string(edges.size()) + ": Number of edges.");
                for (const auto& edge : edges) {
                    auto id = edge->id();
                    auto start = edge->start();
                    auto end = edge->end();
                    auto curve = Curve::create(start, end);
                    auto pair = std::make_pair(id, curve);
                    curves_.insert(pair);
                }
            });

            logger.debug(std::to_string(curves_.size()) + ": Number of curves.");
            logger.debug(std::to_string(data.nodes.size()) + ": Number of nodes.");
            logger.debug(std::to_string(data.S.size()) + ": Number of elements surface matrix.");

            logger.trace("Leaving.");
        }

        static std::vector<el_pot_t>
        reactionPotential_(real_t factor,
                           const std::vector<position_t> &points,
                           const surface_ptr_t& surface) {
            util::Logger logger("simploce::bem_curved_tri_nodes_vertices::reactionPotential_()");
            logger.trace("Entering.");

            logger.debug(std::to_string(points.size()) +
                         ": Number of points for electric potential calculation.");

            logger.trace("Leaving.");
            return std::vector<el_pot_t>{1.0, 0.0};
        }

    }

    bem_calc_ptr_t
    CurvedTriNodesVerticesCalculator::create(const param_ptr_t &param, const surface_ptr_t &surface) {
        return std::make_shared<CurvedTriNodesVerticesCalculator>(param, surface);
    }

    CurvedTriNodesVerticesCalculator::CurvedTriNodesVerticesCalculator(param_ptr_t param, surface_ptr_t surface) :
        param_{std::move(param)}, surface_{std::move(surface)} {
        util::Logger logger("simploce::CurvedTriNodesVerticesCalculator::CurvedTriNodesVerticesCalculator()");
        logger.trace("Entering.");

        if (!param_) {
            util::logAndThrow(logger, "Missing parameters.");
        }
        if (!surface_) {
            util::logAndThrow(logger, "Missing (triangulated) surface.");
        }
        bem_curved_tri_nodes_vertices::initialize_(param_, surface_);

        logger.trace("Leaving.");
    }

    void
    CurvedTriNodesVerticesCalculator::surfaceMatrix() {
        util::Logger logger("simploce::CurvedTriNodesVerticesCalculator::surfaceMatrix()");
        logger.trace("Entering.");

        auto& data = *bem_curved_tri_nodes_vertices::bemData_;

        logger.trace("Leaving.");
    }

    void
    CurvedTriNodesVerticesCalculator::rightHandSide(std::vector<position_t> &positions,
                                                         std::vector<charge_t> &charges) {
        static util::Logger logger("simploce::CurvedTriNodesVerticesCalculator::rightHandSide()");
        logger.trace("Entering.");

        if (positions.size() != charges.size()) {
            util::logAndThrow(logger, "Length charges and positions vector are not the same.");
        }
        if (positions.empty() || charges.empty()) {
            util::logAndThrow(logger, "Positions and/or charges not provided.");
        }

        auto& data = *bem_curved_tri_nodes_vertices::bemData_;
        data.b = rhs::rightHandSide(data.nodes,
                                    positions,
                                    charges,
                                    data.epsSolvent,
                                    data.epsSolute);

        logger.trace("Leaving.");
    }

    void
    CurvedTriNodesVerticesCalculator::solve() {
        static util::Logger logger("simploce::CurvedTriNodesVerticesCalculator::solve()");
        logger.trace("Entering.");
        auto& data = *bem_curved_tri_nodes_vertices::bemData_;
        data.x = data.lu.solve(data.b);
        //std::clog << "Solution: "<< std::endl << data.x << std::endl;
        logger.trace("Leaving.");
    }

    std::vector<el_pot_t>
    CurvedTriNodesVerticesCalculator::reactionPotentialSolute(const std::vector<position_t> &points) {
        util::Logger logger("simploce::CurvedTriNodesVerticesCalculator::reactionPotentialSolute()");
        logger.trace("Entering.");

        logger.debug(std::to_string(points.size()) +
                     ": Number of points for electric potential calculation.");

        auto& data = *bem_curved_tri_nodes_vertices::bemData_;
        auto epsRatio = data.epsRatio;
        auto factor = (epsRatio + 1) / 2.0;
        auto result = bem_curved_tri_nodes_vertices::reactionPotential_(factor, points, surface_);

        logger.trace("Leaving.");
        return std::move(result);
    }

    std::vector<el_pot_t>
    CurvedTriNodesVerticesCalculator::reactionPotentialSolvent(const std::vector<position_t> &points) {
        util::Logger logger("simploce::FlatTriNodesVerticesCalculator::reactionPotentialSolvent()");
        logger.trace("Entering.");

        logger.debug(std::to_string(points.size()) + ": Number of points for electric potential calculation.");
        auto& data = *bem_curved_tri_nodes_vertices::bemData_;
        auto epsRatio = data.epsRatio;
        auto factor = (epsRatio + 1) / (2.0 * epsRatio);
        auto result = bem_curved_tri_nodes_vertices::reactionPotential_(factor, points, surface_);

        logger.trace("Leaving.");
        return std::move(result);
    }

}