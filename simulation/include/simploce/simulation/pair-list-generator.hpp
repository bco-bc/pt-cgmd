/*
 * File:   pair-list-generator.hpp
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
         * @param all All particles.
         * @param free All free particles.
         * @param groups All particle groups.
         * @return Pair lists.
         */
        virtual PairLists
        generate(const std::vector<p_ptr_t>& all,
                 const std::vector<p_ptr_t>& free,
                 const std::vector<pg_ptr_t>& groups) const = 0;
    };
    
}

#endif /* PAIR_LIST_GENERATOR_HPP */

