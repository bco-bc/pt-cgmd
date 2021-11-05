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
            real_t r1 = util::random<real_t>() * MAX;
            real_t r2 = util::random<real_t>() * MAX;
            std::string v1 = util::toString(int(r1));
            std::string v2 = util::toString(int(r2));
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
            return std::nearbyint(val);
        }
    
    }
}