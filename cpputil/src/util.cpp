/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "simploce/util/util.hpp"

namespace simploce {
    namespace util {
    
        int nint(real_t val)
        {
            std::fesetround(FE_TONEAREST);
            return std::nearbyint(val);
        }
    
    }
}