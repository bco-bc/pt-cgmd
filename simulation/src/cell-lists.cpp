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
    
    
    
    // For all particle groups in a single cell.
    template <typename P>
    static typename PairLists<P>::pp_list_cont_t
    forGroupsInOneCell_(const Cell<P>& cell,
                        const bc_ptr_t& bc,
                        const box_ptr_t& box)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        using pp_pair_t = typename PairLists<P>::pp_pair_t;
        
        if ( cell.groups().empty() || cell.groups().empty() ) {
            return pp_list_cont_t{};  // Empty list.
        }
        
        pp_list_cont_t pairList{};

        auto& groups = cell.groups();
        
        for (auto it_i = groups.begin(); it_i != (groups.end() - 1); ++it_i) {
            auto& gi = *it_i;
            const auto particles_i = gi->particles();
            for ( auto it_j = it_i + 1; it_j != groups.end(); ++it_j) {
                auto& gj = *it_j;
                const auto particles_j = gj->particles();
                // Include all particle pairs.
                for (auto pi : particles_i) {
                    for (auto pj : particles_j) {
                        pp_pair_t pair = std::make_pair(pi, pj);
                        pairList.push_back(pair);
                    }
                }
            }
        }
        
        // Done.
        return std::move(pairList);        
    }

    // For all particle groups in two -different- cells.
    template <typename P>
    static typename PairLists<P>::pp_list_cont_t
    forGroupsInTwoCells_(const Cell<P>& cell, 
                         const Cell<P>& neighbor,
                         const bc_ptr_t& bc,
                         const box_ptr_t& box)
    {
        using pp_list_cont_t = typename PairLists<P>::pp_list_cont_t;
        using pp_pair_t = typename PairLists<P>::pp_pair_t;
        
        static real_t rc2 = util::squareCutoffDistance(box);
        
        if ( cell.groups().empty() || cell.groups().empty() ) {
            return pp_list_cont_t{};  // Empty list.
        }

        pp_list_cont_t pairList{};
        
        for (auto& gi: cell.groups()) {
            auto r_gi = gi->position();
            auto& particles_i = gi->particles();
            for (auto gj: neighbor.groups()) {
                auto r_gj = gj->position();                
                auto R = bc->apply(r_gi, r_gj);
                auto R2 = norm2<real_t>(R);
                if ( R2 <= rc2) {
                    // Include all particle pairs.
                    auto& particles_j = gj->particles();
                    for (auto pi :  particles_i) {
                        for (auto pj : particles_j) {
                            pp_pair_t pair = std::make_pair(pi, pj);
                            pairList.push_back(pair);                                
                        }
                    }
                }
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
        
        // Display some information.
        if ( firstTime ) {
            std::clog << "Using cell-based particle pair lists." << std::endl;
            std::clog << "Cell side length: " << grid->sideLength() << std::endl;
            std::clog << "Number of cells: " << grid->numberOfCells() << std::endl;
        }
        
        pp_list_cont_t pairList{};
        
        // Update particle locations in cells.
        grid->place(bc, free, groups);
        
        std::size_t nppff = 0;
        std::size_t nppgg = 0;
        const auto& cells = grid->cells();
        
        // Generate particle pair list. Select particles from cells within 
        // cutoff distance.
        std::vector<std::size_t> numberOfNeighbors{};
        for (auto cell_i = cells.begin(); cell_i != (cells.end() - 1); ++cell_i) {
            auto plCell = forGroupsInOneCell_(*cell_i, bc, box);
            pairList.insert(pairList.end(), plCell.begin(), plCell.end());
            std::size_t nb = 0;
            auto ri = cell_i->position();
            for (auto cell_j = cell_i + 1; cell_j != cells.end(); ++cell_j) {
                auto rj = cell_j->position();
                auto R = bc->apply(ri, rj);
                auto R2 = norm2<real_t>(R);
                if ( R2 <= rc2 ) {
#ifdef _DEBUG
                    std::clog << "Cell pair included: (" 
                              << cell_i->locationAsString() << ") AND ("
                              << cell_j->locationAsString() << ")" << std::endl;
#endif
                    nb += 1;
                    auto plf = forFree_<P>(*cell_i, *cell_j, bc, box);
                    nppff += plf.size();
                    pairList.insert(pairList.end(), plf.begin(), plf.end());
                    auto plg = forGroupsInTwoCells_<P>(*cell_i, *cell_j, bc, box);
                    nppgg += plg.size();
                    pairList.insert(pairList.end(), plg.begin(), plg.end());
                }
            }
            numberOfNeighbors.push_back(nb);
        }
        // Last cell
        auto plCell = forGroupsInOneCell_(*(cells.end() - 1) , bc, box);
        pairList.insert(pairList.end(), plCell.begin(), plCell.end());
        numberOfNeighbors.push_back(0);

        // Display some information.
#ifndef _DEBUG
        if ( firstTime ) {            
#endif
            std::size_t counter = 0;
            std::size_t sum = 0;
            for (auto nb : numberOfNeighbors) {
                sum += nb;
            }
            std::size_t ave = real_t(sum) /  numberOfNeighbors.size();
            std::clog << "Cells:" << std::endl;
            std::clog << "Average number of neighbors: " << ave << std::endl;
            std::size_t ngroups = 0;
            for (auto& cell : cells) {
                
#ifdef _DEBUG
                if ( cell.free().size() > 0 || cell.groups().size() > 0 ) {                    
                    std::clog << "Cell location: (" << cell.locationAsString()
                              << ")," << space;
                    std::clog << "Cell position: " << cell.position() << ", ";
                    std::clog << "Number of free particles: " 
                              << cell.free().size() << ", ";
                    std::clog << "Number of particle groups: " 
                              << cell.groups().size() << ", ";
                    std::clog << "Number of neighbors: "
                              << numberOfNeighbors[counter] << std::endl;                     
                }
#endif
                
                ngroups += cell.groups().size();
                counter += 1;
            }
            assert(ngroups == groups.size());
            std::clog << "Total number of pairs involving free particles: " 
                      << nppff << std::endl;
            std::clog << "Total number of particle-in-group/particle-in-group pairs: "
                      << nppgg << std::endl;
            std::clog << "Total number of particle pairs: "
                      << pairList.size() << std::endl;
            auto total = all.size() * (all.size() - 1) / 2;
            std::clog << "Total number of POSSIBLE particle pairs: "
                      << total << std::endl;
            std::clog << "Fraction (%): " 
                      << real_t(pairList.size()) * 100.0 / total << std::endl;
            firstTime = false;
#ifndef _DEBUG
        }
#endif
        
        // Done.
        return std::move(PairLists<P>(pairList));
    }
    
    CellLists<Atom>::CellLists(const box_ptr_t& box,
                               const bc_ptr_t& bc) :
        box_{box}, bc_{bc}, sideLength_{}
    {  
        sideLength_ = computeSideLength_(box_);
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