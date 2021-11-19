/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 2:09 PM
 */

#ifndef PAIR_LIST_GENERATOR_HPP
#define PAIR_LIST_GENERATOR_HPP

#include "s-types.hpp"
#include "pair-lists.hpp"
#include <vector>

namespace simploce {        
    
    /**
     * Finds all particles pairs in molecular dynamics simulations.
     */
    struct pair_lists_generator {

        virtual ~pair_lists_generator() = default;
        
        /**
         * Generates pair lists.
         * @param particleSystem Particle system.
         * @return Pair lists.
         */
        virtual PairLists
        generate(const p_system_ptr_t& particleSystem) const = 0;
    };
    
}

#endif /* PAIR_LIST_GENERATOR_HPP */

