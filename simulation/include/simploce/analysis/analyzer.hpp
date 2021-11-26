/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 22 September 2019, 15:32
 */

#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include "a-types.hpp"
#include "simploce/particle/p-types.hpp"
#include <vector>

namespace simploce {
    
    /**
     * Performs an analysis on a collection of particles.
     */
    struct Analyzer {
        
        virtual ~Analyzer() = default;
        
        /**
         * Carries out the analysis.
         * @param particleSystem Particle system.
         */
        virtual void perform(const p_system_ptr_t& particleSystem) = 0;
    };
}

#endif /* ANALYZER_HPP */

