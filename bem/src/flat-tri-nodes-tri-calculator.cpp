/*
 * Author: André H. Juffer.
 * Created on 11/07/2022.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/flat-tri-nodes-tri-calculator.hpp"
#include "simploce/bem/bem-data.hpp"
#include "simploce/bem/kernels.hpp"
#include "simploce/bem/rhs.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/face.hpp"
#include "simploce/util/logger.hpp"
#include <memory>

namespace simploce {

    namespace bem_flat_tri_nodes_tri {

        static std::shared_ptr<BEMData> bemData_{};

        /**
         * Assign triangle centers as nodes.
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
                int counter = 0;
                for (auto& face: faces) {
                    BEMData::node_t node;
                    auto pair = face->center();
                    node.r = pair.first;
                    node.nv = pair.second;
                    node.index = counter;
                    node.area = face->area();
                    nodes.emplace_back(node);
                    counter += 1;
                }
                return std::move(nodes);
            });
            return std::move(nodes);
        }

        static void
        initialize_(const param_ptr_t& param, const surface_ptr_t& surface) {
            util::Logger logger("simploce::bem_flat_tri_nodes_tri::initialize_()");
            logger.trace("Entering.");

            bemData_ = std::make_shared<BEMData>(param, surface->numberOfFaces());
            auto& data = *bemData_;

            data.nodes = assignNodes_(surface);
            logger.debug(std::to_string(data.nodes.size()) + ": Number of nodes.");
            logger.debug(std::to_string(surface->numberOfFaces()) + ": Number of triangles.");
            logger.debug(std::to_string(data.S.size()) + ": Number of elements surface matrix.");

            logger.trace("Leaving.");
        }

        static std::vector<el_pot>
        reactionPotential_(real_t factor,
                           const std::vector<position_t> &points,
                           const surface_ptr_t& surface) {
            auto result = surface->doWithAll<std::vector<el_pot>>([factor, points](
                        const std::vector<vertex_ptr_t> &vertices,
                        const std::vector<face_ptr_t> &faces,
                        const std::vector<edge_ptr_t> &edges) {
                std::vector<el_pot> rp(points.size(), 0.0);  // Reaction potentials.
                auto& data = *bem_flat_tri_nodes_tri::bemData_;
                auto epsRatio = data.epsRatio;
                for (auto k = 0; k != points.size(); ++k) {
                    auto& rk = points[k];
                    int counter = 0;
                    for (const auto& face: faces) {
                        auto pair = face->center();
                        auto value = kernels::Lij0(epsRatio, pair.first, pair.second, rk);
                        auto x = data.x[counter];  // Potential at node (triangle face) on boundary.
                        rp[k] += factor * value * face->area()() * x;
                        counter += 1;
                    }
                }
                return std::move(rp);
            });
            return std::move(result);
        }

    }

    bem_calc_ptr_t
    FlatTriNodesTriCalculator::create(const param_ptr_t &param, const surface_ptr_t &surface) {
        return std::make_shared<FlatTriNodesTriCalculator>(param, surface);
    }

    FlatTriNodesTriCalculator::FlatTriNodesTriCalculator(param_ptr_t param, surface_ptr_t surface) :
    param_{std::move(param)}, surface_{std::move(surface)} {
        util::Logger logger{"simploce::FlatTriNodesTriCalculator::FlatTriNodesTriCalculator()"};
        logger.trace("Entering.");

        if (!param_) {
            util::logAndThrow(logger, "Missing parameters.");
        }
        if (!surface_) {
            util::logAndThrow(logger, "Missing (triangulated) surface.");
        }
        bem_flat_tri_nodes_tri::initialize_(param_, surface_);

        logger.trace("Leaving.");
    }

    void
    FlatTriNodesTriCalculator::surfaceMatrix() {
        util::Logger logger("simploce::FlatTriNodesTriCalculator::surfaceMatrix()");
        logger.trace("Entering");

        surface_->doWithAll<void>([&logger] (const std::vector<vertex_ptr_t>& vertices,
                                                  const std::vector<face_ptr_t>& faces,
                                                  const std::vector<edge_ptr_t>& edges) {
            auto& data = *bem_flat_tri_nodes_tri::bemData_;

            auto& S = data.S;
            auto size = S.size();
            logger.debug(std::to_string(size) + ": Size (number of elements) of surface matrix.");

            // Compute elements of S.
            auto& nodes = data.nodes;
            logger.debug(std::to_string(nodes.size()) + ": Number of nodes.");
            S.setZero();
            for (auto i = 0; i != nodes.size(); ++i) {
                S(i,i) = 1.0;
            }

            auto epsRatio = data.epsRatio;
            for (auto i = 0; i != nodes.size(); ++i) {
                auto& node_i = nodes[i];
                auto index_i = node_i.index;
                auto r0 = node_i.r;
                for (auto j = 0; j != nodes.size(); ++j) {
                    auto& node_j = nodes[j];
                    auto value = i != j ?
                                        kernels::Lij0(epsRatio, node_j.r, node_j.nv, r0) :
                                        0.0;
                    value *= node_j.area();
                    S(index_i,node_j.index) -= value;
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
    FlatTriNodesTriCalculator::rightHandSide(std::vector<position_t> &positions, std::vector<charge_t> &charges) {
        util::Logger logger("simploce::FlatTriNodesTriCalculator::rightHandSide");
        logger.trace("Entering");

        if (positions.size() != charges.size()) {
            util::logAndThrow(logger, "Length charges and positions vector are not the same.");
        }
        if (positions.empty() || charges.empty()) {
            util::logAndThrow(logger, "Positions and/or charges not provided.");
        }

        auto& data = *bem_flat_tri_nodes_tri::bemData_;
        data.b = rhs::rightHandSide(data.nodes,positions, charges, data.epsSolvent, data.epsSolute);

        logger.trace("Leaving");
    }

    void
    FlatTriNodesTriCalculator::solve() {
        static util::Logger logger("simploce::FlatTriNodesTriCalculator::solve()");
        logger.trace("Entering.");
        auto& data = *bem_flat_tri_nodes_tri::bemData_;
        data.x = data.lu.solve(data.b);
        //std::clog << "Solution: "<< std::endl << data.x << std::endl;
        logger.trace("Leaving");
    }

    std::vector<el_pot>
    FlatTriNodesTriCalculator::reactionPotentialSolute(const std::vector<position_t> &points) {
        static util::Logger logger("simploce::FlatTriNodesTriCalculator::reactionPotentialSolute()");
        logger.trace("Entering.");

        auto data = *bem_flat_tri_nodes_tri::bemData_;
        auto epsRatio = data.epsRatio;
        auto factor = (epsRatio + 1) / 2.0;
        auto result = bem_flat_tri_nodes_tri::reactionPotential_(factor, points, surface_);

        logger.trace("Leaving.");
        return std::move(result);
    }

    std::vector<el_pot>
    FlatTriNodesTriCalculator::reactionPotentialSolvent(const std::vector<position_t> &points) {
        static util::Logger logger("simploce::FlatTriNodesTriCalculator::reactionPotentialSolvent()");
        logger.trace("Entering.");

        auto data = *bem_flat_tri_nodes_tri::bemData_;
        auto epsRatio = data.epsRatio;
        auto factor = (epsRatio + 1) / (2.0 * epsRatio);
        auto result = bem_flat_tri_nodes_tri::reactionPotential_(factor, points, surface_);

        logger.trace("Leaving.");
        return std::move(result);
    }
}
