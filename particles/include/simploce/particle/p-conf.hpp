/*
 * File:   p-conf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 12, 2019, 4:12 PM
 */

#ifndef P_CONF_HPP
#define P_CONF_HPP

#include "p-types.hpp"
#include "simploce/conf/u-conf.hpp"

namespace simploce {
    namespace conf {
    
        /**
         * Width of a (particle, spec) name output field in an output stream.
         */
        const int NAME_WIDTH = 10;
        
        /**
         * No changes allowed in particle's protonation state.
         */
        const std::size_t NOT_PROTONATABLE = 0;

        /**
         * Changes allowed in particle's protonation state.
         */
        const std::size_t PROTONATABLE = 1;

    }
}

#endif /* P_CONF_HPP */

