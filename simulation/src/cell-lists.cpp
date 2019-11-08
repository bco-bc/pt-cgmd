/*
 * The MIT License
 *
 * Copyright 2019 juffer.
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

#include "simploce/simulation/cell-lists.hpp"
#include "simploce/simulation/grid.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/simulation/sim-util.hpp"
#include "simploce/particle/particle-group.hpp"
#include <memory>
#include <vector>
#include <thread>
#include <cassert>

namespace simploce {
    
    /**
     * Returns grid element side length.
     * @param box Box.
     * @return Side length.
     */
    static length_t 
    sideLength_(const box_ptr_t& box)
    {
        return util::cutoffDistance(box);
    }
    
    
    
    // Between free and other free particles, and free and particles in groups.
    template <typename P>
    static typename PairLists<P>::pp_list_cont_t
    forFree_(const std::vector<std::shared_ptr<P>>& free, 
             const Cell<P>& cell,
             const bc_ptr_t& bc,
             real_t rc2)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        
        pp_list_cont_t pairList{};
        
        /*
        for (auto f : free) {
            auto index = f->index();
            position_t rf = f->position();
            for (const auto& p : cell.free()) {
                // Avoid double counting, or interacting with itself.
                if ( p->index() > index ) {
                    position_t rp = p->position();
                    auto dv = bc->apply(rf,rp);
                    auto r2 = norm2<real_t>(dv);
                    if ( r2 <= rc2 ) {
                        // Include this pair.
                        p_pair_t pair = std::make_pair(f, p);
                        pairList.push_back(pair);
                    }
                }
            }
            for (const auto& pg: cell.groups()) {
                auto rg = pg->position();
                auto dv = bc->apply(rf, rg);
                auto r2 = norm2<real_t>(dv);
                if ( r2 <= rc2) {
                    // Include all interactions between free particle and the 
                    // particles in the given group.
                    for ( auto p : pg->particles()) {
                        p_pair_t pair = std::make_pair(f, p);
                        pairList.push_back(pair);                        
                    }
                }
            }
        }
        */
        
        return std::move(pairList);
    }

    // Between particles in groups.
    template <typename P>
    static typename PairLists<P>::pp_list_cont_t
    forGroups_(const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups,
               const Cell<P>& cell,
               const bc_ptr_t& bc,
               real_t rc2)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        
        pp_list_cont_t pairList{};
        
        /*
        for (auto gi: groups) {
            position_t ri = gi->position();
            const auto& particles_i = gi->particles();
            for (auto gj: cell.groups()) {
                // Avoid particles in the same group.
                if (gj != gi ) {
                   position_t rj = gj->position();
                    auto dv = bc->apply(ri, rj);
                    auto r2 = norm2<real_t>(dv);
                    if ( r2 <= rc2) {
                        // Include all interactions between particles of this 
                        // group pair.
                        const auto& particles_j = gj->particles();
                        for (auto pi :  particles_i) {
                            auto index = pi->index();
                            for (auto pj : particles_j) {
                                // Avoid double counting.
                                if ( pj->index() > index ) {
                                    p_pair_t pair = std::make_pair(pi, pj);
                                    pairList.push_back(pair);
                                }
                            }
                        }
                    }
                }
            }
        }
         */
        
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
        using grid_t = Grid<P>;
        using grid_ptr_t = typename grid_t::grid_ptr_t;
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
            
        static bool firstTime = true;
        static auto sideLength = sideLength_(box);        
        static grid_ptr_t grid = grid_t::make(box, sideLength);
        static real_t rc2 = util::squareCutoffDistance(box);
        
        pp_list_cont_t pairList{};
        
        grid->clear();
        grid->place(bc, free, groups);
        
        const auto& cells = grid->cells();
        for (auto& cell : cells) {
            const auto& free = cell.free();
            const auto& groups = cell.groups();
            auto neighbors = grid->neighbors(cell.location());
            for (auto& neighbor : neighbors) {
                auto pl = forFree_<P>(free, neighbor, bc, rc2);
                pairList.insert(pairList.end(), pl.begin(), pl.end());
                pl = forGroups_<P>(groups, neighbor, bc, rc2);
                pairList.insert(pairList.end(), pl.begin(), pl.end());
            }
        }
        
        if ( firstTime ) {
            std::clog << "Total number of particle pairs: "
                      << pairList.size() << std::endl;
            std::clog << "Total number of POSSIBLE particle pairs: "
                      << all.size() * (all.size() - 1) / 2 << std::endl;
            firstTime = false;
        }
        
        // Done.
        return std::move(PairLists<P>(pairList));
    }
    
    CellLists<Atom>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
    
    PairLists<Atom>
    CellLists<Atom>::generate(const std::vector<atom_ptr_t>& all,
                              const std::vector<atom_ptr_t>& free,
                              const std::vector<atom_group_ptr_t>& groups) const    
    {
        return std::move(makePairLists_<Atom>(box_, bc_, all, free, groups));
    }
    
    CellLists<Bead>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
    
    PairLists<Bead> 
    CellLists<Bead>::generate(const std::vector<bead_ptr_t>& all,
                              const std::vector<bead_ptr_t>& free,
                              const std::vector<bead_group_ptr_t>& groups) const    {
        return std::move(makePairLists_<Bead>(box_, bc_, all, free, groups));
    }    
    
}