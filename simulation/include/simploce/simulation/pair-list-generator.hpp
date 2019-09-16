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
 * File:   pair-list-generator.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 2:09 PM
 */

#ifndef PAIR_LIST_GENERATOR_HPP
#define PAIR_LIST_GENERATOR_HPP

#include "stypes.hpp"
#include "bc.hpp"
#include "simploce/particle/particle-group.hpp"
#include <vector>
#include <utility>
#include <thread>
#include <memory>

namespace simploce {
    
    /**
     * @param P Particle type.
     * @param BC Boundary condition type.
     */
    template <typename P>
    class ParticlePairListGenerator {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle pair type.
         */
        using p_pair_t = std::pair<p_ptr_t, p_ptr_t>;
        
        /**
         * Particle pair lists type (single list).
         */
        using p_pair_list_t = std::vector<p_pair_t>;
        
        /**
         * Particle pair lists type (multiple lists).
         */
        using p_pair_lists_t = std::vector<p_pair_list_t>;
        
        /**
         * Particle group type.
         */
        using pg_t = ParticleGroup<P>;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<pg_t>;
        
        /**
         * Constructor.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @param dt Time step.
         */        
        ParticlePairListGenerator(const box_ptr_t& box,
                                  const bc_ptr_t& bc) : box_{box}, bc_{bc} {}
        
        /**
         * Generates particle pairs lists.
         * @return List, multiple or just one. Depends on the number of available threads.
         */
        std::vector<p_pair_list_t> generate(const std::vector<p_ptr_t>& all,
                                            const std::vector<p_ptr_t>& free,
                                            const std::vector<pg_ptr_t>& groups) const;
        
        /**
         * Generates particle pairs lists.
         * @param All particles.
         * @param free All free particles.
         * @param groups Particle groups.
         * @return List, multiple or just one. Depends on the number of available threads.
         */
        p_pair_lists_t operator () (const std::vector<p_ptr_t>& all,
                                    const std::vector<p_ptr_t>& free,
                                    const std::vector<pg_ptr_t>& groups) const;
                
        
    private:
        
        const length_t RINCLUDE_{2.5};
        
        // For any particle collections.
        p_pair_list_t forParticles_(const std::vector<p_ptr_t>& particles) const;
        
        // For particles within groups.
        p_pair_list_t forGroups_(const std::vector<pg_ptr_t>& groups) const;
        
        // For free particles and particles in groups.
        p_pair_list_t forFreeAndGroups_(const std::vector<p_ptr_t>& free,
                                        const std::vector<pg_ptr_t>& groups) const;
        
        /**
         * Returns square of suitable cutoff distance.
         */
        real_t rc2_() const;
        
        box_ptr_t box_;
        bc_ptr_t bc_;
    };
    
    
    template <typename P>    
    typename ParticlePairListGenerator<P>::p_pair_lists_t
    ParticlePairListGenerator<P>::generate(const std::vector<p_ptr_t>& all,
                                           const std::vector<p_ptr_t>& free,
                                           const std::vector<pg_ptr_t>& groups) const
    {
        // Multiple particle pair lists. Its size will be the number of available 
        // threads.
        p_pair_lists_t pairlists{};
    
        // Prepare new particle pair list.
        p_pair_list_t pairlist{};        // All pairs.
        
        p_pair_list_t plFree = this->forParticles_(free);
        pairlist.insert(pairlist.end(), plFree.begin(), plFree.end());
        
        p_pair_list_t plGroups = this->forGroups_(groups);
        pairlist.insert(pairlist.end(), plGroups.begin(), plGroups.end());
       
        p_pair_list_t plFreeGroups = this->forFreeAndGroups_(free, groups);
        pairlist.insert(pairlist.end(), plFreeGroups.begin(), plFreeGroups.end());
        
        // Prepare sublists.                       
        static const std::size_t nlists = std::thread::hardware_concurrency();        
        std::size_t npairsSubList = pairlist.size() / nlists;                                                              
        std::size_t counter = 0;      
        for (std::size_t k = 0; k != nlists; ++k) {
            p_pair_list_t single{};  // One list of particle pairs.
            std::size_t n = 0;
            while (counter != pairlist.size() && n != npairsSubList) {
                single.push_back(pairlist[counter]);
                counter += 1;
                n += 1;
            }
            pairlists.push_back(single);          // Store list.
        }
    
        // Add remaining particles pairs, if any, to the last set.
        if ( counter < pairlist.size() ) {
            p_pair_list_t& last =  *(pairlists.end() - 1);
            while ( counter != pairlist.size() ) {
                last.push_back(pairlist[counter]);
                counter += 1;
            }
        }
        
        // Done.
        return std::move(pairlists);
    }
    
    template <typename P>
    typename ParticlePairListGenerator<P>::p_pair_lists_t 
    ParticlePairListGenerator<P>::operator () (const std::vector<p_ptr_t>& all,
                                               const std::vector<p_ptr_t>& free,
                                               const std::vector<pg_ptr_t>& groups) const
    {
        return std::move(this->generate(all, free, groups));
    }
    
    template <typename P>
    real_t ParticlePairListGenerator<P>::rc2_() const
    {
        length_t boxSize = box_->size();
        length_t rc = (RINCLUDE_() > boxSize() ? boxSize : RINCLUDE_);
        real_t rc2 = rc * rc;
        return rc2;        
    }
    
    template <typename P> 
    typename ParticlePairListGenerator<P>::p_pair_list_t 
    ParticlePairListGenerator<P>::forParticles_(const std::vector<p_ptr_t>& particles) const
    {
        if ( particles.empty() ) {
            return p_pair_list_t{};  // Empty list.
        }
        
        real_t rc2 = this->rc2_();
        p_pair_list_t pairs{};
        
        // For all pairs in the set.
        for (auto iter_i = particles.begin(); iter_i != particles.end() - 1; ++iter_i) {
            p_ptr_t pi = *iter_i;
            position_t ri = pi->position();
            for (auto iter_j = iter_i + 1; iter_j != particles.end(); ++iter_j) {
                p_ptr_t pj = *iter_j;
                position_t rj = pj->position();
                dist_vect_t rij = bc_->apply(ri, rj);
                real_t R2 = norm2<real_t>(rij);
                if ( R2 <= rc2 ) {
                    // Include this pair.
                    p_pair_t pair = std::make_pair(pi, pj);
                    pairs.push_back(pair);
                }
            }
        }

        // Done.
        return pairs;
    }
    
    template <typename P> 
    typename ParticlePairListGenerator<P>::p_pair_list_t 
    ParticlePairListGenerator<P>::forGroups_(const std::vector<pg_ptr_t>& groups) const
    {
        using p_ptr_cont_t = typename pg_t::p_ptr_cont_t;
        
        if ( groups.empty() ) {
            return p_pair_list_t{};  // Empty list.
        }
            
        real_t rc2 = this->rc2_();
        p_pair_list_t pairs{};

        for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
            pg_t& pg_i = **iter_i;
            position_t ri = pg_i.position();
            p_ptr_cont_t particles_i = pg_i.particles();
            for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                pg_t& pg_j = **iter_j;
                position_t rj = pg_j.position();
                dist_vect_t rij = bc_->apply(ri, rj);
                real_t R2 = norm2<real_t>(rij);
                if ( R2 <= rc2 ) {
                    // Include all interactions between particles for this pair.
                    p_ptr_cont_t particles_j = pg_j.particles();
                    for (auto pi :  particles_i) {
                        for (auto pj : particles_j) {
                            p_pair_t pair = std::make_pair(pi, pj);
                            pairs.push_back(pair);
                        }
                    }
                }
            }
        }

        // Done.
        return pairs;        
    }
    
    template <typename P> 
    typename ParticlePairListGenerator<P>::p_pair_list_t 
    ParticlePairListGenerator<P>::forFreeAndGroups_(const std::vector<p_ptr_t>& free,
                                                    const std::vector<pg_ptr_t>& groups) const
    {
        using p_ptr_cont_t = typename pg_t::p_ptr_cont_t;

        if ( free.empty() || groups.empty() ) {
            return p_pair_list_t{};  // Empty list.
        }
        
        real_t rc2 = this->rc2_();
        p_pair_list_t pairs{};
        
        for (auto piter = free.begin(); piter != free.end(); ++piter) {
            p_ptr_t pi = *piter;
            position_t r = pi->position();
            for (auto giter = groups.begin(); giter != groups.end(); ++giter) {
                pg_t& pg = **giter;
                position_t rg = pg.position();
                dist_vect_t rpg = bc_->apply(r, rg);
                real_t R2 = norm2<real_t>(rpg);
                if ( R2 <= rc2 ) {
                    // Include all interactions between free particle and the 
                    // particles in the given group.
                    p_ptr_cont_t particles = pg.particles();
                    for (auto pj : particles) {
                        p_pair_t pair = std::make_pair(pi, pj);
                        pairs.push_back(pair);
                    }
                }
            }
        }
    
        // Done.
        return pairs;
    }
}

#endif /* PAIR_LIST_GENERATOR_HPP */

