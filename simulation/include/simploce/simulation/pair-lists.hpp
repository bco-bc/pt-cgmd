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
     * Holds all particle/particle, particle/group and group/group pair lists. Here,
     * a particle always refers to a -free- particle.
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
         * Particle/particle group pair type.
         */
        using pg_pair_t = std::pair<p_ptr_t, pg_ptr_t>;
        
        /**
         * Particle group/particle group pair type.
         */
        using gg_pair_t = std::pair<pg_ptr_t, pg_ptr_t>;
                
        /**
         * Particle/particle pair list container type.
         */
        using pp_list_cont_t = std::vector<pp_pair_t>;
                
        /**
         * Particle/particle group pair list container type.
         */
        using pg_list_cont_t = std::vector<pg_pair_t>;
        
        /**
         * Particle group/particle group pair list container type.
         */
        using gg_list_cont_t = std::vector<gg_pair_t>;
        
        /**
         * Default constructor. Empty lists.
         */
        PairLists();
        
        /**
         * Constructor.
         * @param ppPairList Particle/particle pair list.
         * @param pgpairList Particle/particle group pair list.
         * @param ggPairList Particle group/particle group pair list.
         */
        PairLists(const pp_list_cont_t& ppPairList,
                  const pg_list_cont_t& pgPairList,
                  const gg_list_cont_t& ggPairList);
        
        /**
         * Returns particle/particle pair list.
         * @return pair list.
         */
        const pp_list_cont_t& particlePairList() const { 
            return ppPairList_; 
        }
        
        /**
         * Returns particle/particle group pair list.
         * @return Pair list.
         */
        const pg_list_cont_t& particleGroupPairList() const { 
            return pgPairList_; 
        }
        
        /**
         * Returns particle group/particle group pair list.
         * @return pair list.
         */
        const gg_list_cont_t& groupPairList() const { 
            return ggPairList_; 
        }
        
        /**
         * Were the pair lists altered.
         * @return Result.
         */
        bool isModified() const { return modified_; }
        
    private:
        
        template <typename PP>
        friend class Interactor;
        
        /**
         * Signals whether the pair lists were modified.
         * @param modified If true, lists were updated, otherwise not.
         */
        void updated_(bool modified) { modified_ = modified; }
        
        pp_list_cont_t ppPairList_;
        pg_list_cont_t pgPairList_;
        gg_list_cont_t ggPairList_;
        bool modified_;
    };
    
    template <typename P>
    PairLists<P>::PairLists() :
        ppPairList_{}, pgPairList_{}, ggPairList_{}, modified_{true}
    {
    }
        
    template <typename P>
    PairLists<P>::PairLists(const pp_list_cont_t& ppPairList,
                            const pg_list_cont_t& pgPairList,
                            const gg_list_cont_t& ggPairList) :
        ppPairList_{ppPairList}, pgPairList_{pgPairList}, 
        ggPairList_{ggPairList}, modified_{true}
    {
    }

}

#endif /* PAIR_LISTS_HPP */

