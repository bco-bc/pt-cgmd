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
    
    }
}