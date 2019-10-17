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
    
    static real_t rcutoff2_(const box_ptr_t& box)
    {
        static length_t rc = util::cutoffDistance(box);
        static real_t rc2 = rc * rc;
        
        return rc2;
    }
    
    /**
     * For any collection of particles.
     * @param particles Particles.
     * @return Pair list.
     */
    template <typename P> 
    typename DistanceLists<P>::p_pair_list_t 
    forParticles_(const box_ptr_t& box,
                  const bc_ptr_t& bc,
                  const std::vector<typename DistanceLists<P>::p_ptr_t>& particles)
    {
        using p_pair_list_t = typename DistanceLists<P>::p_pair_list_t;
        
        static real_t rc2 = rcutoff2_(box);
        
        if ( particles.empty() ) {
            return p_pair_list_t{};  // Empty list.
        }
        
        // Pair list.
        p_pair_list_t pairs{};
        
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
                    pairs.push_back(pair);
                }
            }
        }

        // Done.
        return pairs;
    }
    
    template <typename P> 
    typename DistanceLists<P>::p_pair_list_t 
    forGroups_(const box_ptr_t& box,
               const bc_ptr_t& bc,
               const std::vector<typename DistanceLists<P>::pg_ptr_t>& groups)
    {        
        using p_pair_list_t = typename DistanceLists<P>::p_pair_list_t;
        
        static real_t rc2 =  rcutoff2_(box);

        if ( groups.empty() ) {
            return p_pair_list_t{};  // Empty list.
        }
            
        // Pair list.
        p_pair_list_t pairs{};

        // For all pairs of groups.
        for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
            const auto& pg_i = **iter_i;
            position_t ri = pg_i.position();
            const auto& particles_i = pg_i.particles();
            for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                const auto& pg_j = **iter_j;
                position_t rj = pg_j.position();
                dist_vect_t R = bc->apply(ri, rj);
                real_t R2 = norm2<real_t>(R);
                if ( R2 < rc2 ) {
                    // Include all interactions between particles for this pair
                    // of groups.
                    const auto& particles_j = pg_j.particles();
                    for (auto pi :  particles_i) {
                        for (auto pj : particles_j) {
                            auto pair = std::make_pair(pi, pj);
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
    typename DistanceLists<P>::p_pair_list_t 
    forParticlesAndGroups_(const box_ptr_t& box,
                           const bc_ptr_t& bc,
                           const std::vector<typename DistanceLists<P>::p_ptr_t>& particles,
                           const std::vector<typename DistanceLists<P>::pg_ptr_t>& groups)
    {
        using p_pair_list_t = typename DistanceLists<P>::p_pair_list_t;
        
        static real_t rc2 =  rcutoff2_(box);

        if ( particles.empty() || groups.empty() ) {
            return p_pair_list_t{};  // Empty list.
        }

        // Pair list.
        p_pair_list_t pairs{};
        
        for (auto piter = particles.begin(); piter != particles.end(); ++piter) {
            const auto& pi = *piter;
            position_t r = pi->position();
            for (auto giter = groups.begin(); giter != groups.end(); ++giter) {
                const auto& pg = **giter;
                position_t rg = pg.position();
                dist_vect_t R = bc->apply(r, rg);
                real_t R2 = norm2<real_t>(R);
                if ( R2 < rc2 ) {
                    // Include all interactions between the free particle and the 
                    // particles in the given group.
                    const auto& particles = pg.particles();
                    for (auto pj : particles) {
                        auto pair = std::make_pair(pi, pj);
                        pairs.push_back(pair);
                    }
                }
            }
        }
    
        // Done.
        return pairs;
    }
        
    template <typename P>
    static std::vector<typename DistanceLists<P>::p_pair_list_t>
    makePairLists_(const box_ptr_t& box,
                   const bc_ptr_t& bc,
                   const std::vector<typename DistanceLists<P>::p_ptr_t>& all,
                   const std::vector<typename DistanceLists<P>::p_ptr_t>& free,
                   const std::vector<typename DistanceLists<P>::pg_ptr_t>& groups)
    {
        using p_pair_list_t = typename DistanceLists<P>::p_pair_list_t;
        
        static bool firstTime = true;
        if ( firstTime ) {
            std::clog << "Using particles pair lists based on distances between particles." 
                      << std::endl;
            std::clog << "Cutoff distance: " << util::cutoffDistance(box) << std::endl;
            firstTime = false;
        }
        
        // Prepare new particle pair list.
        // All pairs.
        p_pair_list_t pairList{};
        
        auto plFree = forParticles_<P>(box, bc, free);
        pairList.insert(pairList.end(), plFree.begin(), plFree.end());
        
        auto plGroups = forGroups_<P>(box, bc, groups);
        pairList.insert(pairList.end(), plGroups.begin(), plGroups.end());
        
        auto pl = forParticlesAndGroups_<P>(box, bc, free, groups);
        pairList.insert(pairList.end(), pl.begin(), pl.end());
        
        // Prepare sublists, the number of which will depend on the number of 
        // hardware threads available.               
        std::vector<p_pair_list_t> pairLists{};
        static const std::size_t nlists = std::thread::hardware_concurrency();        
        std::size_t npairsSubList = pairList.size() / nlists;                                                              
        std::size_t counter = 0;      
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
            p_pair_list_t& last =  *(pairLists.end() - 1);
            while ( counter != pairList.size() ) {
                last.push_back(pairList[counter]);
                counter += 1;
            }
        }
        
        // Done.
        return pairLists;
    }
    
    DistanceLists<Atom>::DistanceLists(const box_ptr_t& box,
                                       const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
        
    std::vector<typename DistanceLists<Atom>::p_pair_list_t> 
    DistanceLists<Atom>::generate(const std::vector<p_ptr_t>& all,
                                  const std::vector<p_ptr_t>& free,
                                  const std::vector<pg_ptr_t>& groups) const 
    {
        return std::move(makePairLists_<Atom>(box_, bc_, all, free, groups));
    }
    
    DistanceLists<Bead>::DistanceLists(const box_ptr_t& box,
                                       const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
        
    std::vector<typename DistanceLists<Bead>::p_pair_list_t> 
    DistanceLists<Bead>::generate(const std::vector<p_ptr_t>& all,
                                  const std::vector<p_ptr_t>& free,
                                  const std::vector<pg_ptr_t>& groups) const 
    {
        return makePairLists_<Bead>(box_, bc_, all, free, groups);
    }
    
}