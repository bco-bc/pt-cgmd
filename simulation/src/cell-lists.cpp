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
    computeSideLength_(const box_ptr_t& box)
    {
        return 0.5 * util::cutoffDistance(box);
        return util::cutoffDistance(box);
    }
    
    
    
    // Between free and neighboring free particles, and free and particles in 
    // neighboring groups.
    template <typename P>
    static typename PairLists<P>::pp_list_cont_t
    forFree_(const Cell<P>& cell, 
             const Cell<P>& neighbor,
             const bc_ptr_t& bc,
             const box_ptr_t& box)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        using pp_pair_t = typename PairLists<P>::pp_pair_t;
        
        static real_t rc2 = util::squareCutoffDistance(box);
        
        pp_list_cont_t pairList{};
        
        for (auto& fi : cell.free()) {
            auto index_i = fi->index();
            position_t ri = fi->position();
            for (auto& fj : neighbor.free()) {                
                auto index_j = fj->index();
                if ( index_j > index_i ) {
                    position_t rj = fj->position();
                    auto R = bc->apply(ri,rj);
                    auto R2 = norm2<real_t>(R);
                    if ( R2 <= rc2 ) {
                        // Include this pair.
                        pp_pair_t pair = std::make_pair(fi, fj);
                        pairList.push_back(pair);
                    }
                }
            }
            for (const auto& pg: neighbor.groups()) {
                auto rg = pg->position();
                auto R = bc->apply(ri, rg);
                auto R2 = norm2<real_t>(R);
                if ( R2 <= rc2) {
                    // Include all pairs between given free particle and particles
                    // in group.
                    for ( auto pj : pg->particles()) {
                        auto index_j = pj->index();
                        if ( index_j > index_i ) { 
                            pp_pair_t pair = std::make_pair(fi, pj);
                            pairList.push_back(pair);                        
                       }
                    }
                }
            }
        }
        
        return std::move(pairList);
    }

    // Between particles in different groups.
    template <typename P>
    static typename PairLists<P>::pp_list_cont_t
    forGroups_(const Cell<P>& cell, 
               const Cell<P>& neighbor,
               const bc_ptr_t& bc,
               const box_ptr_t& box)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        using pp_pair_t = typename PairLists<P>::pp_pair_t;
        
        //static real_t rc2 = util::squareCutoffDistance(box);

        pp_list_cont_t pairList{};
        
        for (auto& gi: cell.groups()) {
            //position_t ri = gi->position();
            auto& particles_i = gi->particles();
            for (auto gj: neighbor.groups()) {
                //position_t rj = gj->position();
                /*
                auto R = bc->apply(ri, rj);
                auto R2 = norm2<real_t>(R);
                if ( R2 <= rc2) {
                    // Include all interactions between the particles of this 
                    // group pair.
                 */
                    auto& particles_j = gj->particles();
                    for (auto pi :  particles_i) {
                        //auto index_i = pi->index();
                        for (auto pj : particles_j) {
                            //auto index_j = pj->index();
                            // Avoid double counting.
                            //if ( index_j > index_i ) { 
                                pp_pair_t pair = std::make_pair(pi, pj);
                                pairList.push_back(pair);                                
                            //}
                        }
                    }
                //}
            }
        }
        
        return std::move(pairList);        
    }
    
    template <typename P>
    static PairLists<P>
    makePairLists_(const box_ptr_t& box,
                   const bc_ptr_t& bc,
                   const length_t& sideLength,
                   const std::vector<std::shared_ptr<P>>& all,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {
        using grid_t = Grid<P>;
        using grid_ptr_t = typename grid_t::grid_ptr_t;
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
            
        static bool firstTime = true;                
        static grid_ptr_t grid = grid_t::make(box, sideLength);
        static auto rc2 = util::squareCutoffDistance(box);
        
        if ( firstTime ) {
            std::clog << "Using cell-based particle pair lists." << std::endl;
            std::clog << "Cell side length: " << grid->sideLength() << std::endl;
            std::clog << "Number of cells: " << grid->numberOfCells() << std::endl;
        }
        
        pp_list_cont_t pairList{};
        
        // Update particle location in cells.
        grid->clear();
        grid->place(bc, free, groups);
        
        // Generate particle pair list.
        const auto& cells = grid->cells();
        for (auto it_i = cells.begin(); it_i != cells.end() - 1; ++it_i) {
            auto& cell_i = *it_i;
            auto ri = cell_i.position();
            for (auto it_j = it_i + 1; it_j != cells.end(); ++it_j) {
                auto& cell_j = *it_j;
                auto rj = cell_j.position();
                auto R = bc->apply(ri, rj);
                auto R2 = norm2<real_t>(R);
                if ( R2 <= rc2 ) {
                    auto plf = forFree_<P>(cell_i, cell_j, bc, box);
                    pairList.insert(pairList.end(), plf.begin(), plf.end());
                    auto plg = forGroups_<P>(cell_i, cell_j, bc, box);
                    pairList.insert(pairList.end(), plg.begin(), plg.end());
                }
            }
        }
        /*
        for (auto& cell : cells) {
            auto neighbors = grid->neighbors(cell.location());
            for (auto& neighbor : neighbors) {
                auto pl = forFree_<P>(cell, neighbor, bc, box);
                pairList.insert(pairList.end(), pl.begin(), pl.end());
                pl = forGroups_<P>(cell, neighbor, bc, box);
                pairList.insert(pairList.end(), pl.begin(), pl.end());
            }
        }
        */
        if ( firstTime ) {
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
    
    CellLists<Atom>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}, sideLength_{computeSideLength_(box)}
    {        
    }
    
    PairLists<Atom>
    CellLists<Atom>::generate(const std::vector<atom_ptr_t>& all,
                              const std::vector<atom_ptr_t>& free,
                              const std::vector<atom_group_ptr_t>& groups) const    
    {
        return std::move(makePairLists_<Atom>(box_, 
                                              bc_, 
                                              sideLength_, 
                                              all, 
                                              free, 
                                              groups));
    }
    
    void 
    CellLists<Atom>::sideLength(const length_t& sideLength)
    {
        sideLength_ = sideLength;
    }
    
    CellLists<Bead>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}, sideLength_{computeSideLength_(box)}
    {        
    }
    
    PairLists<Bead> 
    CellLists<Bead>::generate(const std::vector<bead_ptr_t>& all,
                              const std::vector<bead_ptr_t>& free,
                              const std::vector<bead_group_ptr_t>& groups) const    {
        return std::move(makePairLists_<Bead>(box_, 
                                              bc_, 
                                              sideLength_, 
                                              all, 
                                              free, 
                                              groups));
    }    
    
    void 
    CellLists<Bead>::sideLength(const length_t& sideLength)
    {
        sideLength_ = sideLength;
    }
    
}