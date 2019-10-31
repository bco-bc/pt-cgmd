/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   pair-lists.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 21, 2019, 3:43 PM
 */

#ifndef PAIR_LISTS_HPP
#define PAIR_LISTS_HPP

#include "simploce/particle/particle-group.hpp"
#include <utility>
#include <memory>
#include <vector>

namespace simploce {
    
    /**
     * Holds all particle pairs for non-bonded interactions.
     * @param P Particle type.
     */
    template <typename P>
    class PairLists {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<ParticleGroup<P>>;
        
        /**
         * Particle/particle pair type.
         */
        using pp_pair_t = std::pair<p_ptr_t, p_ptr_t>;
        
        /**
         * Particle/particle pair list container type.
         */
        using pp_list_cont_t = std::vector<pp_pair_t>;
                
        /**
         * Default constructor. Empty list.
         */
        PairLists();
        
        /**
         * Constructor.
         * @param ppPairList Particle/particle pair list.
         */
        PairLists(const pp_list_cont_t& pairList);
        
        /**
         * Returns particle/particle pair list.
         * @return pair list.
         */
        const pp_list_cont_t& particlePairList() const { 
            return pairList_; 
        }
               
        /**
         * Was the pair list altered.
         * @return Result.
         */
        bool isModified() const { return modified_; }
        
    private:
        
        template <typename PP>
        friend class Interactor;
        
        /**
         * Signals whether the pair list was modified.
         * @param modified If true, list was updated, otherwise not.
         */
        void updated_(bool modified);
        
        pp_list_cont_t pairList_;
        bool modified_;
    };
    
    template <typename P>
    PairLists<P>::PairLists() :
        pairList_{}, modified_{true}
    {
    }
        
    template <typename P>
    PairLists<P>::PairLists(const pp_list_cont_t& pairList) :
        pairList_{pairList}, modified_{true}
    {
    }
        
    template <typename P>
    void 
    PairLists<P>::updated_(bool modified)
    {
        modified_ = modified;
    }

}

#endif /* PAIR_LISTS_HPP */

