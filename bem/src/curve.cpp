/*
 * Author: André H. Juffer.
 * Created on 15/07/2022
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/curve.hpp"
#include "simploce/surface/vertex.hpp"
#include <tuple>
#include <array>

namespace simploce {
    namespace curve {

        /**
         * Returns constant vectors a, b, c, and d of curve between two vertices.
         * @param start Start vertex.
         * @param end End vertex.
         * @return Constant vectors a, b, c, d.
         */
        static std::tuple<position_t, position_t, position_t, position_t>
        calculateCurve_(const vertex_ptr_t& start, const vertex_ptr_t& end) {
            auto positionStart = start->position();
            auto positionEnd = end->position();
            auto R = positionEnd - positionStart;
            auto normalStart = start->normal();
            auto normalEnd = end->normal();
            auto imp1 = inner<real_t>(R, normalStart);
            auto imp2=inner<real_t>(R, normalEnd);
            auto imp3 = inner<real_t>(normalEnd, normalEnd);
            auto imp4 = inner<real_t>(normalStart, normalStart);
            auto imp5 = inner<real_t>(normalStart, normalEnd);

            // Calculate constant vectors.
            auto f1 = (imp1 * imp3 + 0.5 * imp2 * imp5) / (imp4 * imp3 / 3.0 - (imp5 * imp5) / 12.0);
            auto f2= (-1.0 * imp2 - f1 * imp5 / 6.0) / (imp3 / 3.0);
            position_t a, b, c, d;
            for (auto k = 0; k != 3; ++k) {
                a[k] = positionStart[k];
                c[k] = 0.5 * f1 * normalStart[k];
                d[k] = (f2 * normalEnd[k] - f1 * normalStart[k]) / 6.0;
                b[k] = positionEnd[k] - a[k] - c[k] - d[k];
            }
            return std::make_tuple(a, b, c, d);
        }

    }

    curve_ptr_t
    Curve::create(const vertex_ptr_t &start, const vertex_ptr_t &end) {
        return std::make_shared<Curve>(start, end);
    }

    Curve::Curve(vertex_ptr_t start, vertex_ptr_t end) :
        start_{std::move(start)}, end_{std::move(end)} {
        auto constantVectors = curve::calculateCurve_(start_, end_);
        a_ = std::get<0>(constantVectors);
        b_ = std::get<1>(constantVectors);
        c_ = std::get<2>(constantVectors);
        d_ = std::get<3>(constantVectors);
    }

    std::pair<position_t, normal_t>
    Curve::point(real_t t) {
        auto tt= t * t;
        auto ttt = t * tt;

        // Position.
        position_t r;
        for (auto k = 0; k != 3; ++k) {
            r[k] = a_[k] + b_[k] * t + c_[k] * tt + d_[k] * ttt;
        }

        std::array<real_t,3> der1{0.0, 0.0, 0.0}, der2{0.0, 0.0, 0.0};
        real_t imp{0}, l_der1{0};
        for (auto k = 0; k != 3; ++k) {
            der1[k] = b_[k] + 2.0 * c_[k] * t + 3.0 * d_[k] * tt;
            der2[k] = 2.0 * c_[k] + 6.0 * d_[k] * t;
            l_der1 += der1[k] * der1[k];
            imp += der1[k] * der2[k];
        }

        // Curvature vector (cv) and square of length curvature vector (l_cv_2).
        std::array<real_t, 3> cv{0.0, 0.0, 0.0};
        real_t l_cv_2{0};
        for (auto k = 0; k != 3; ++k) {
            cv[k] = (der2[k] - imp * der1[k] / l_der1) / l_der1;
            l_cv_2 += cv[k] * cv[k];
        }

        // Unit normal vector.
        normal_t normal;
        auto ave = (start_->normal() + end_->normal()) / 2.0;
        if (l_cv_2 > 0.0) {
            // Curve has curvature.
            auto l_cv = std::sqrt(l_cv_2);  // Length of curvature vector.
            for (auto k = 0; k != 3; ++k) {
                normal[k] = cv[k] / l_cv;
            }
        } else {
            // Take the average of start and end vertex normal vector.
            normal = normal / norm<real_t>(ave);
        }

        // Ensure same orientation as start and end vertex normal.
        normal = inner<real_t>(normal, ave) > 0.0 ? normal : -1.0 * normal;

        return std::move(std::make_pair(r, normal));
    }
}