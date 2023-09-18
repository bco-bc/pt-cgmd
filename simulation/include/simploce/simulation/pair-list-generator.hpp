/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 2:09 PM
 */

#ifndef PAIR_LIST_GENERATOR_HPP
#define PAIR_LIST_GENERATOR_HPP

#include "simploce/types/s-types.hpp"
#include "pair-list.hpp"
#include <vector>

namespace simploce {        
    
    /**
     * Finds non-bonded particles pairs in particle systems.
     */
    struct pair_list_generator {

        /**
         * Finds non-bonded particle pairs.
         * @param particleSystem Particle system.
         * @return Particle pairs.
         */
        virtual std::vector<PairList::p_pair_t>
        generate(const p_system_ptr_t& particleSystem) const = 0;

    };
    
}

#endif /* PAIR_LIST_GENERATOR_HPP */

