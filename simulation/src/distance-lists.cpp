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
 * File:   simple-lists.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 12:21 PM
 */

#include "simploce/simulation/distance-lists.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/sim-util.hpp"
#include <vector>
#include <utility>
#include <thread>

namespace simploce {
    
    static const std::size_t MIN_NUMBER_OF_PARTICLES = 100;        
    
    static real_t rcutoff2_(const box_ptr_t& box)
    {
        static length_t rc = util::cutoffDistance(box);
        static real_t rc2 = rc * rc;
        
        return rc2;
    }
    
    /**
     * For -any- collection of particles.
     */
    template <typename P> 
    typename PairLists<P>::pp_list_cont_t
    forParticles_(const box_ptr_t& box,
                  const bc_ptr_t& bc,
                  const std::vector<std::shared_ptr<P>>& particles)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        
        static real_t rc2 = rcutoff2_(box);
        
        if ( particles.empty() ) {
            return pp_list_cont_t{};  // Empty pair list.
        }
        
        // Particle/particle pair list.
        pp_list_cont_t ppPairList{};
        
        // For all pairs in the set.
        for (auto iter_i = particles.begin(); iter_i != particles.end() - 1; ++iter_i) {
            const auto& pi = *iter_i;
            position_t ri = pi->position();
            for (auto iter_j = iter_i + 1; iter_j != particles.end(); ++iter_j) {
                const auto& pj = *iter_j;
                position_t rj = pj->position();
                dist_vect_t R = bc->apply(ri, rj);
                real_t R2 = norm2<real_t>(R);                
                if ( R2 < rc2 ) {
                    // Include this pair.
                    auto pair = std::make_pair(pi, pj);
                    ppPairList.push_back(pair);
                }
            }
        }

        // Done.
        return std::move(ppPairList);
    }
    
    template <typename P> 
    typename PairLists<P>::gg_list_cont_t
    forGroups_(const box_ptr_t& box,
               const bc_ptr_t& bc,
               const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {    
        using gg_list_cont_t = typename PairLists<P>::gg_list_cont_t;
        using gg_pair_t = typename PairLists<P>::gg_pair_t;
        
        static real_t rc2 =  rcutoff2_(box);

        if ( groups.empty() ) {
            return gg_list_cont_t{};  // Empty list.
        }
            
        // Group/Group pair list.
        gg_list_cont_t ggPairList{};

        // For all pairs of different groups.
        for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
            const auto pg_i = *iter_i;
            position_t ri = pg_i->position();            
            for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                const auto pg_j = *iter_j;
                position_t rj = pg_j->position();
                dist_vect_t R = bc->apply(ri, rj);
                real_t R2 = norm2<real_t>(R);
                if ( R2 < rc2 ) {
                    gg_pair_t pair = std::make_pair(pg_i, pg_j);
                    ggPairList.push_back(pair);
                }
            }
        }

        // Done.
        return std::move(ggPairList);        
    }
    
    template <typename P> 
    typename PairLists<P>::pg_list_cont_t
    forParticlesAndGroups_(const box_ptr_t& box,
                           const bc_ptr_t& bc,
                           const std::vector<std::shared_ptr<P>>& particles,
                           const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {
        using pg_list_cont_t = typename PairLists<P>::pg_list_cont_t;
        using pg_pair_t = typename PairLists<P>::pg_pair_t;
        
        static real_t rc2 =  rcutoff2_(box);

        if ( particles.empty() || groups.empty() ) {
            return pg_list_cont_t{};  // Empty list.
        }

        // Pair list.
        pg_list_cont_t pgPairList{};
        
        for (auto p : particles) {
            position_t rp = p->position();
            for (auto g : groups) {
                position_t rg = g->position();
                dist_vect_t R = bc->apply(rp, rg);
                real_t R2 = norm2<real_t>(R);
                if ( R2 < rc2 ) {
                    pg_pair_t pair = std::make_pair(p, g);
                    pgPairList.push_back(pair);
                }
            }
        }
    
        // Done.
        return std::move(pgPairList);
    }
        
    template <typename P>
    static PairLists<P>
    makePairLists_(const box_ptr_t& box,
                   const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& all,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {
        static bool firstTime = true;
        if ( firstTime ) {
            std::clog << "Using particles pair lists based on distances between particles." 
                      << std::endl;
            std::clog << "Cutoff distance: " << util::cutoffDistance(box) << std::endl;
            std::clog << "Total number of particles: " << all.size() << std::endl;
            std::clog << "Total number of free particles: " << free.size() << std::endl;
            std::clog << "Total number of particle groups: " << groups.size() << std::endl;
        }
        
        // Prepare new particle pair list.        
        auto ppPairList = forParticles_<P>(box, bc, free);       
        auto pgPairList = forParticlesAndGroups_<P>(box, bc, free, groups);
        auto ggPairList = forGroups_<P>(box, bc, groups);        
        
        /*
        // Prepare sublists, the number of which will depend on the number of 
        // hardware threads available.               
        std::vector<p_pair_list_t> pairLists{};
        if ( all.size() > conf::MIN_NUMBER_OF_PARTICLES ) {
            std::size_t counter = 0;      
            static const std::size_t nlists = std::thread::hardware_concurrency();        
            std::size_t npairsSubList = pairList.size() / nlists;                                                              
            for (std::size_t k = 0; k != nlists; ++k) {
                p_pair_list_t single{};  // One sublist of particle pairs.
                std::size_t n = 0;
                while (counter != pairList.size() && n != npairsSubList) {
                    single.push_back(pairList[counter]);
                    counter += 1;
                    n += 1;
                }
                pairLists.push_back(single);          // Store list.
            }
    
            // Add remaining particles pairs, if any, to the last set.
            if ( counter < pairList.size() ) {
                p_pair_list_t& last = *(pairLists.end() - 1);
                while ( counter != pairList.size() ) {
                    last.push_back(pairList[counter]);
                    counter += 1;
                }
            }
        } else {
            pairLists.push_back(pairList);
        }
        */
        
        if ( firstTime ) {
            std::clog << "Number of free particle pairs: " << ppPairList.size() 
                      << std::endl;
            std::clog << "Number of free particle/particle group pairs: "
                      << pgPairList.size() << std::endl;
            std::clog << "Number of particle group pairs: "
                      << ggPairList.size() << std::endl;
            firstTime = false;
        }

        // Done.
        return std::move(PairLists<P>(ppPairList, pgPairList, ggPairList));
    }
    
    DistanceLists<Atom>::DistanceLists(const box_ptr_t& box,
                                       const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
        
    PairLists<Atom>
    DistanceLists<Atom>::generate(const std::vector<atom_ptr_t>& all,
                                  const std::vector<atom_ptr_t>& free,
                                  const std::vector<atom_group_ptr_t>& groups) const 
    {
        return makePairLists_<Atom>(box_, bc_, all, free, groups);
    }
    
    DistanceLists<Bead>::DistanceLists(const box_ptr_t& box,
                                       const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
        
    PairLists<Bead>
    DistanceLists<Bead>::generate(const std::vector<bead_ptr_t>& all,
                                  const std::vector<bead_ptr_t>& free,
                                  const std::vector<bead_group_ptr_t>& groups) const 
    {
        return makePairLists_<Bead>(box_, bc_, all, free, groups);
    }
    
}