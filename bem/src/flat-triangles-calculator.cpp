/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 20:58.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/flat-triangles-calculator.hpp"
#include "simploce/bem/bem-data.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/vertex.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include "simploce/units/units-mu.hpp"
#include <vector>
#include <stdexcept>
#include <memory>

namespace simploce {

    std::shared_ptr<BEMData> bemData_{};

    bem_calc_ptr_t
    FlatTrianglesCalculator::create(const param_ptr_t &param, const surface_ptr_t &surface) {
        return std::make_shared<FlatTrianglesCalculator>(param, surface);
    }

    FlatTrianglesCalculator::FlatTrianglesCalculator(param_ptr_t param, surface_ptr_t surface) :
            param_(std::move(param)), surface_{std::move(surface)} {
        util::Logger logger("simploce::FlatTrianglesCalculator::FlatTrianglesCalculator()");
        logger.trace("Entering");

        if (!param_) {
            throw std::domain_error("FlatTrianglesCalculator: Missing parameters.");
        }
        if (!surface_) {
            throw std::domain_error("FlatTrianglesCalculator: Missing (triangulated) surface.");
        }
        auto ka = param_->get<real_t>("bem.solvent.ka");
        auto numberOfFaces = surface_->numberOfFaces();
        auto nCol= ka > 0.0 ? 2 * numberOfFaces : numberOfFaces;
        bemData_ = std::make_shared<BEMData>(param_, nCol);

        logger.debug(std::to_string(ka) + ": Value of inverse Debye length.");
        logger.debug(std::to_string(numberOfFaces) + ": Number of triangles.");
        auto square = numberOfFaces * numberOfFaces;
        logger.debug(std::to_string(square) + ": Square number of triangles.");

        logger.trace("Leaving.");
    }

    void
    FlatTrianglesCalculator::surfaceMatrix() {
        util::Logger logger("simploce::FlatTrianglesCalculator::surfaceMatrix()");
        logger.trace("Entering");

        surface_->doWithAll<void>([this, &logger] (const std::vector<vertex_ptr_t>& vertices,
                                                        const std::vector<face_ptr_t>& faces,
                                                        const std::vector<edge_ptr_t>& edges) {

            auto& data = *bemData_;
            auto nCol = data.S.size();
            logger.debug(std::to_string(nCol) + ": Dimension of surface matrix.");

            // Compute elements of S.
            auto surface = this->surface_;
            data.S.setZero();

            // LU-decomposition.
            data.lu.compute(data.S);
        });

        // Done.
        logger.trace("Leaving.");
    }

    void
    FlatTrianglesCalculator::rightHandSide(std::vector<position_t>& positions, std::vector<charge_t>& charges) {
        util::Logger logger("simploce::FlatTrianglesCalculator::rightHandSide");
        logger.trace("Entering");

        if (positions.size() != charges.size()) {
            util::logAndThrow(logger, "Length charges and positions vector are not the same.");
        }
        if (positions.empty() || charges.empty()) {
            util::logAndThrow(logger, "Positions and/or charges not provided.");
        }

        // Right-hand-side of Sx=b.
        surface_->doWithAll<void >([&positions, &charges, &logger] (const std::vector<vertex_ptr_t>& vertices,
                                                                         const std::vector<face_ptr_t>& faces,
                                                                         const std::vector<edge_ptr_t>& edges) {
            auto& b = bemData_->b;
            b.setZero();
            auto FOUR_PI_E0 = units::mu<real_t>::FOUR_PI_E0;
            auto t1 = 2.0 / ((1.0 + bemData_->epsRatio) * FOUR_PI_E0 * bemData_->epsSolute);
            auto N = BEMData::index_t(vertices.size());
            logger.debug(std::to_string(N) + ": Number of collocation points.");
            auto ka = bemData_->ka;
            if (ka > 0.0) {
                logger.debug(std::to_string(bemData_->ka) + ": Inverse Debye length. Non-zero ionic strength.");
                auto t2 = t1 * bemData_->epsRatio;
                for (auto i = 0; i != charges.size(); ++i) {
                    auto Q = charges[i]();
                    auto r = positions[i];

                    for (auto j = 0; j != N; ++j) {
                        auto& vertex = vertices[j];
                        auto r0 = vertex->position();
                        auto unv0 = vertex->normal();
                        auto R = r - r0;
                        auto dis = norm<real_t>(R);
                        auto dis3 = dis * dis * dis;
                        auto imp = inner<real_t>(R, unv0);
                        b[j] += t1 * Q / dis;
                        b[j + N] += t2 * Q * imp / dis3;
                    }
                }
            } else {
                logger.debug(std::to_string(bemData_->ka) + ": Inverse Debye length. Zero ionic strength.");
                for (auto i = 0; i != charges.size(); ++i) {
                    auto Q = charges[i]();
                    auto r = positions[i];
                    for (auto j = 0; j != N; ++j) {
                        auto& vertex = vertices[j];
                        auto r0 = vertex->position();
                        auto R = r - r0;
                        auto dis = norm<real_t>(R);
                        b[j] += t1 * Q / dis;
                    }
                }
            }
        });

        logger.trace("Leaving");
    }

    void
    FlatTrianglesCalculator::solve() {
        auto& data = *bemData_;
        data.x = data.lu.solve(data.b);
    }

    std::vector<el_pot>
    FlatTrianglesCalculator::electricPotentials(const std::vector<position_t> &points) {
        util::Logger logger("simploce::FlatTrianglesCalculator::electricPotentials()");
        logger.trace("Entering.");

        auto result =
                surface_->doWithAll<std::vector<el_pot>>([&logger, points](
                        const std::vector<vertex_ptr_t> &vertices,
                        const std::vector<face_ptr_t> &faces,
                        const std::vector<edge_ptr_t> &edges) {
               logger.debug(std::to_string(points.size()) + ": Number of points for electric potential calculation.");
               auto& data = *bemData_;
               std::vector<el_pot> ep(points.size(), 0.0);
               return std::move(ep);
           });

        logger.trace("Leaving.");
        return std::move(result);
    }

}