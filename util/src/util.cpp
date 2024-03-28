/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "simploce/util/util.hpp"
#include "boost/uuid/uuid_generators.hpp"
#include "boost/lexical_cast.hpp"

namespace simploce {
    namespace util {

        // @see https://en.wikipedia.org/wiki/Gaussian_quadrature
        using GaussianQuadrature = std::pair<std::vector<real_t>, std::vector<real_t>>;

        static GaussianQuadrature
        makeEightPointGaussianQuadrature() {
            std::vector<real_t> points(8);
            std::vector<real_t> weights(8);

            points[0] = -0.960289856497536;
            points[1] = -0.796666477413627;
            points[2] = -0.525532409916329;
            points[3] = -0.183434642495650;
            points[4] = 0.960289856497536;
            points[5] = 0.796666477413627;
            points[6] = 0.525532409916329;
            points[7] = 0.183434642495650;

            weights[0] = 0.101228536200376;
            weights[1] = 0.222381034453374;
            weights[2] = 0.313706645877887;
            weights[3] = 0.362683783378362;
            weights[4] = 0.101228536200376;
            weights[5] = 0.222381034453374;
            weights[6] = 0.313706645877887;
            weights[7] = 0.362683783378362;

            return std::make_pair(points, weights);
        }

        static GaussianQuadrature QUADRATURE = makeEightPointGaussianQuadrature();

        /*
        id_t generateId() {
            return boost::uuids::random_generator()();
        }

        id_t generateId() {
            static int ID = 0;
            return ++ID;
        }
         */

        id_t generateId() {
            static real_t MAX = 99999;
            real_t r1 = util::random() * MAX;
            real_t r2 = util::random() * MAX;
            std::string v1 = std::to_string(int(r1));
            std::string v2 = std::to_string(int(r2));
            return v1 + "-" + v2;
        }

        /*
        std::string toString(const simploce::id_t& id) {
            return boost::uuids::to_string(id);
        }
         */

        id_t toId(const std::string& v) {
            return boost::lexical_cast<id_t>(v);
        }
    
        int nint(real_t val)
        {
            std::fesetround(FE_TONEAREST);
            return int(std::nearbyint(val));
        }

        real_t randomUniform(real_t min, real_t max) {
            std::random_device rd;
            std::mt19937_64 gen(rd());
            std::uniform_real_distribution<real_t> dis(min, max);
            auto value = dis(gen);
            return value;
        }

        real_t random() {
            return randomUniform(0.0,1.0);
        }

        dist_vect_t randomUnit() {
            auto x = util::random();
            x = util::random() > 0.5 ? x : -x;
            auto y = util::random();
            y = util::random() > 0.5 ? y : -y;
            auto z = util::random();
            z = util::random() > 0.5 ? z : -z;
            dist_vect_t uv{x, y, z};
            auto length = simploce::norm<real_t>(uv);
            return uv / length;
        }

        /**
         * Uses a 8 point Gaussian quadrature rule.
         *
        real_t
        integrate(real_t FUNCTION(const real_t x),
                  real_t a,
                  real_t b) {
            static auto points = QUADRATURE.first;
            static auto weights = QUADRATURE.second;

            auto t1 = (b - a) / 2.0;
            auto t2 = (a + b) / 2.0;
            auto dx_dp = (b - a) / 2.0;

            real_t result = 0.0;
            for (size_t k = 0; k != points.size(); ++k) {
                auto& p = points[k];
                auto& w = weights[k];
                auto x = t1 * p + t2;
                result += FUNCTION(x) * w;
            }
            result *= dx_dp;

            return result;
        }
         */

        /**
         * Uses a 8 point Gaussian quadrature rule.
         */
        real_t
        integrate(const std::function<real_t (real_t x)>& function,
                  real_t a,
                  real_t b) {
            static auto points = QUADRATURE.first;
            static auto weights = QUADRATURE.second;

            auto t1 = (b - a) / 2.0;
            auto t2 = (a + b) / 2.0;
            auto dx_dp = (b - a) / 2.0;

            real_t result = 0.0;
            for (size_t k = 0; k != points.size(); ++k) {
                auto& p = points[k];
                auto& w = weights[k];
                auto x = t1 * p + t2;
                result += function(x) * w;
            }
            result *= dx_dp;

            return result;
        }
    
    }
}