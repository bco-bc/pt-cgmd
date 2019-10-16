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
#include "simploce/particle/particle-group.hpp"
#include <memory>
#include <vector>
#include <thread>
#include <cassert>

namespace simploce {
    
    // Default cutoff distance for non bonded interactions.
    static length_t RCUTOFF_DISTANCE_{2.5};  // nm.
    
    /**
     * Returns cutoff distance.
     * @param box Box.
     * @return Cutoff Distance, always <= 0.5 * box size.
     */
    static length_t 
    rc_(const box_ptr_t& box)
    {
        length_t halve = 0.5 * box->size();
        length_t rc = RCUTOFF_DISTANCE_() > halve() ? halve : RCUTOFF_DISTANCE_;
        return rc;
    }
    
    /**
     * Returns square of cutoff distance.
     * @param box Box.
     * @return Square.
     */
    static real_t 
    rc2_(const box_ptr_t& box)
    {
        length_t rc = rc_(box);
        return rc * rc;
    }
                
    // Between free and other free particles, and free and particles in groups.
    template <typename P>
    static typename CellLists<P>::p_pair_list_t
    forFree_(const std::vector<typename CellLists<P>::p_ptr_t>& free, 
             const Cell<P>& cell,
             const bc_ptr_t& bc,
             real_t rc2)
    {
        using p_pair_list_t = typename CellLists<P>::p_pair_list_t;
        using p_pair_t = typename CellLists<P>::p_pair_t;
        using p_ptr_cont_t = typename ParticleGroup<P>::p_ptr_cont_t;
        
        p_pair_list_t pairList{};
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
                    p_ptr_cont_t particles = pg->particles();
                    for ( auto p : particles) {
                        // Avoid double counting, or interacting with itself..
                        if ( p->index() > index ) {
                            p_pair_t pair = std::make_pair(f, p);
                            pairList.push_back(pair);
                        }
                    }
                }
            }
        }
        
        return std::move(pairList);
    }

    // Between particles in groups.
    template <typename P>
    static typename CellLists<P>::p_pair_list_t
    forGroups_(const std::vector<typename CellLists<P>::pg_ptr_t>& groups,
               const Cell<P>& cell,
               const bc_ptr_t& bc,
               real_t rc2)
    {
        using p_pair_list_t = typename CellLists<P>::p_pair_list_t;
        using p_ptr_cont_t = typename ParticleGroup<P>::p_ptr_cont_t;
        using p_pair_t = typename CellLists<P>::p_pair_t;
        
        p_pair_list_t pairList{};
        for (auto gi: groups) {
            position_t ri = gi->position();
            auto particles_i = gi->particles();
            for (const auto& gj: cell.groups()) {
                // Avoid particles in the same group.
                if (gj != gi ) {
                   position_t rj = gj->position();
                    auto dv = bc->apply(ri, rj);
                    auto r2 = norm2<real_t>(dv);
                    if ( r2 <= rc2) {
                        // Include all interactions between particles of this 
                        // group pair.
                        p_ptr_cont_t particles_j = gj->particles();
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
        
        return std::move(pairList);        
    }
    
    template <typename P>
    static std::vector<typename CellLists<P>::p_pair_list_t>
    makePairLists_(const box_ptr_t& box,
                   const bc_ptr_t& bc,
                   const std::vector<typename CellLists<P>::p_ptr_t>& all,
                   const std::vector<typename CellLists<P>::p_ptr_t>& free,
                   const std::vector<typename CellLists<P>::pg_ptr_t>& groups)
    {
        using cell_t = Cell<P>;
        using grid_t = Grid<P>;
        using p_pair_list_t = typename CellLists<P>::p_pair_list_t;
            
        static grid_t grid{};
        static real_t rc2;
        
        // Set up if still required.
        static bool setup = false;
        if ( !setup ) {
            auto rc = rc_(box);
            grid = Grid<P>::make(box, rc);
            rc2 = rc2_(box);
            setup = true;
        }
        
        // Determine full particle pair list.
        p_pair_list_t full{};
        grid.clear();
        grid.place(bc, free, groups);
        const auto& cells = grid.cells();
        for (auto iter = cells.begin(); iter != cells.end(); ++iter) {
            const cell_t& cell = *iter;
            const auto& free = cell.free();
            const auto& groups = cell.groups();
            auto neighbors = grid.neighbors(cell);
            for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
                cell_t& neighbor = *it;
                p_pair_list_t pl = forFree_<P>(free, neighbor, bc, rc2);
                full.insert(full.end(), pl.begin(), pl.end());
                pl = forGroups_<P>(groups, neighbor, bc, rc2);
                full.insert(full.end(), pl.begin(), pl.end());
            }
        }
        
        // Prepare sublists, the number of which depends on number of 
        // hardware threads available.
        std::vector<p_pair_list_t> pairLists{};   // Sublists.
        static const std::size_t nlists = std::thread::hardware_concurrency();        
        std::size_t npairsSubList = full.size() / nlists;                                                              
        std::size_t counter = 0;      
        for (std::size_t k = 0; k != nlists; ++k) {
            p_pair_list_t single{};  // One list of particle pairs.
            std::size_t n = 0;
            while (counter != full.size() && n != npairsSubList) {
                single.push_back(full[counter]);
                counter += 1;
                n += 1;
            }
            pairLists.push_back(single);          // Store list.
        }
    
        // Add remaining particles pairs, if any, to the last sublist.
        if ( counter < full.size() ) {
            p_pair_list_t& last =  *(pairLists.end() - 1);
            while ( counter != full.size() ) {
                last.push_back(full[counter]);
                counter += 1;
            }
        }
        
        // Done.
        return std::move(pairLists);        
    }
    
    CellLists<Atom>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
    
    std::vector<CellLists<Atom>::p_pair_list_t> 
    CellLists<Atom>::generate(const std::vector<p_ptr_t>& all,
                              const std::vector<p_ptr_t>& free,
                              const std::vector<pg_ptr_t>& groups) const
    {
        return std::move(makePairLists_<Atom>(box_, bc_, all, free, groups));
    }
    
    CellLists<Bead>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}
    {        
    }
    
    std::vector<CellLists<Bead>::p_pair_list_t> 
    CellLists<Bead>::generate(const std::vector<p_ptr_t>& all,
                              const std::vector<p_ptr_t>& free,
                              const std::vector<pg_ptr_t>& groups) const
    {
        return std::move(makePairLists_<Bead>(box_, bc_, all, free, groups));
    }
    
    
}