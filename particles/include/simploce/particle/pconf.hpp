/*
 * File:   pconf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 12, 2019, 4:12 PM
 */

#ifndef PCONF_HPP
#define PCONF_HPP

#include "ptypes.hpp"
#include "simploce/util/uconf.hpp"

namespace simploce {
    namespace conf {
    
        /**
         * Width of a (particle, spec) name output field in an output stream.
         */
        const int NAME_WIDTH = 10;
        
        /**
         * Mass of a proton (u)
         */
        const mass_t MASS_PROTON = 1.007276466879;
        
        /**
         * Charge of a proton (e).
         */
        const charge_t CHARGE_PROTON = 1.0;
        
        /**
         * No changes in protonation state.
         */
        const std::size_t NOT_PROTONATABLE = 0;
        
        /**
         * Discrete changes in particle's protonation state.
         */
        const std::size_t DISCRETELY_PROTONATABLE = 1;
        
        /**
         * Continouos changes in particle's protonation state.
         */
        const std::size_t CONTINUOUSLY_PROTONATABLE = 2;
        
    }
}

#endif /* CONF_HPP */

