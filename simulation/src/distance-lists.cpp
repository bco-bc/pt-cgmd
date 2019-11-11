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
        
        static real_t rc2 = util::squareCutoffDistance(box);
        
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
                if ( R2 <= rc2 ) {
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
    typename PairLists<P>::pp_list_cont_t
    forGroups_(const box_ptr_t& box,
               const bc_ptr_t& bc,
               const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {    
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        using pp_pair_t = typename PairLists<P>::pp_pair_t;
        
        static real_t rc2 =  util::squareCutoffDistance(box);

        if ( groups.empty() ) {
            return pp_list_cont_t{};  // Empty list.
        }
            
        // Pair list.
        pp_list_cont_t pairList{};

        // For all particle group pairs.
        for (auto iter_i = groups.begin(); iter_i != groups.end() - 1; ++iter_i) {
            const auto gi = *iter_i;
            auto r_gi = gi->position();
            const auto particles_i = gi->particles();
            for (auto iter_j = iter_i + 1; iter_j != groups.end(); ++iter_j) {
                auto gj = *iter_j;
                auto r_gj = gj->position();
                auto R = bc->apply(r_gi, r_gj);
                auto R2 = norm2<real_t>(R);
                if ( R2 <= rc2 ) {
                    // Include all particle pairs.
                    const auto particles_j = gj->particles();
                    for (auto pi : particles_i) {
                        for (auto pj : particles_j) {
                            pp_pair_t pair = std::make_pair(pi, pj);
                            pairList.push_back(pair);
                        }
                    }
                }
            }
        }
        
        // Done.
        return std::move(pairList);        
    }
    
    template <typename P> 
    typename PairLists<P>::pp_list_cont_t
    forParticlesAndGroups_(const box_ptr_t& box,
                           const bc_ptr_t& bc,
                           const std::vector<std::shared_ptr<P>>& particles,
                           const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        using pp_pair_t = typename PairLists<P>::pp_pair_t;
        
        static real_t rc2 =  util::squareCutoffDistance(box);

        if ( particles.empty() || groups.empty() ) {
            return pp_list_cont_t{};  // Empty list.
        }

        // Pair list.
        pp_list_cont_t pairList{};
        
        for (auto pi : particles) {
            position_t ri = pi->position();
            for (auto g : groups) {
                if ( !g->contains(pi) ) {
                    for (auto pj : g->particles()) {
                        auto rj = pj->position();
                        auto R = bc->apply(ri, rj);
                        real_t R2 = norm2<real_t>(R);
                        if ( R2 <= rc2 ) {
                            pp_pair_t pair = std::make_pair(pi, pj);
                            pairList.push_back(pair);
                        }
                    }
                }
            }
        }
    
        // Done.
        return std::move(pairList);
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
            std::clog << "Using distance-based particle pair lists." 
                      << std::endl;
            std::clog << "Cutoff distance: " << util::cutoffDistance(box) << std::endl;
        }
        
        // Prepare new particle pair list.        
        auto pairList = forParticles_<P>(box, bc, free);
        auto ppSize = pairList.size();
        auto fgPairList = forParticlesAndGroups_<P>(box, bc, free, groups);
        auto fgSize = fgPairList.size();
        pairList.insert(pairList.end(), fgPairList.begin(), fgPairList.end());
        auto ggPairList = forGroups_<P>(box, bc, groups);
        auto ggSize = ggPairList.size();
        pairList.insert(pairList.end(), ggPairList.begin(), ggPairList.end());
        
        if ( firstTime ) {
            std::clog << "Number of free-particle/free-particle pairs: " 
                      << ppSize << std::endl;
            std::clog << "Number of free-particle/particle-in-group pairs: "
                      << fgSize << std::endl;
            std::clog << "Number of particle-in-group/particle-in-group pairs: "
                      << ggSize << std::endl;
            std::clog << "Total number of particle pairs: "
                      << pairList.size() << std::endl;
            auto total = all.size() * (all.size() - 1) / 2;
            std::clog << "Total number of POSSIBLE particle pairs: "
                      << total << std::endl;
            std::clog << "Fraction (%): " 
                      << real_t(pairList.size()) * 100.0 / total << std::endl;
            firstTime = false;
        }

        // Done.
        return std::move(PairLists<P>(pairList));
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
        return std::move(makePairLists_<Atom>(box_, bc_, all, free, groups));
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
        return std::move(makePairLists_<Bead>(box_, bc_, all, free, groups));
    }
    
}