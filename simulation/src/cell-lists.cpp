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
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle-group.hpp"
#include <boost/multi_array.hpp>
#include <tuple>
#include <memory>
#include <vector>
#include <thread>

namespace simploce {
    
    // Cell.
    template <typename P>
    struct Cell {
        using p_ptr_t = typename CellLists<P>::p_ptr_t;
        using pg_ptr_t = typename CellLists<P>::pg_ptr_t;
        using indices_t = std::tuple<std::size_t, std::size_t, std::size_t>;
        
        Cell() : r_{}, all_{}, free_{}, groups_{}, indices_{} {}
        
        Cell(const position_t& r, 
             const indices_t& indices) : 
            r_{r}, all_{}, free_{}, groups_{}, indices_{indices} {}
        
        position_t r_;
        std::vector<p_ptr_t> all_;
        std::vector<p_ptr_t> free_;
        std::vector<pg_ptr_t> groups_;
        indices_t indices_;
    };
    
    // Cutoff distance for non bonded interactions.
    static length_t RCUTOFF_DISTANCE_{2.5};  // nm.
    
    /**
     * Returns cutoff distance.
     * @param box Box.
     * @return Distance. 
     */
    static length_t 
    rc_(const box_ptr_t& box)
    {
        length_t rc = 0.5 * box->size();
        return RCUTOFF_DISTANCE_() > rc() ? rc : RCUTOFF_DISTANCE_;        
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
    
    template <typename P>
    static std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, std::size_t, std::size_t>
    neighbors_(const Cell<P>& cell, std::size_t size)
    {
        const auto& indices = cell.indices_;
        auto l = std::get<0>(indices);
        auto i_i = l > 0 ? l - 1 : l;
        auto i_f = l < size - 1 ? size : l;
        l = std::get<1>(indices);
        auto j_i = l > 0 ? l - 1 : l;
        auto j_f = l < size - 1 ? size : l; 
        l = std::get<2>(indices);
        auto k_i = l > 0 ? l - 1 : l;
        auto k_f = l < size - 1 ? size : l; 
        return std::make_tuple(i_i, i_f, j_i, j_f, k_i, k_f);
    }
    
    /**
     * Makes cells
     * @param box Simulation box.
     * @param rc Cutoff distance.
     * @return Cells, cell edge length, and size per dimension.
     */
    template <typename P>
    static std::tuple<boost::multi_array<Cell<P>, 3>, length_t, std::size_t>
    makeCells_(const box_ptr_t& box,
               const length_t& rc)
    {
        using cell_t = Cell<P>;
        using cells_t = boost::multi_array<cell_t, 3>;
        using indices_t = typename Cell<P>::indices_t;
        
        auto edgeLength = 0.5 * rc;
        std::size_t n = box->edgeLength() / edgeLength();
        cells_t cells(boost::extents[n][n][n]);
        for (std::size_t i = 0; i != n; ++i) {
            real_t x = (i + 0.5) * edgeLength();
            for (std::size_t j = 0; j != n; ++j) {
                real_t y = (j + 0.5) * edgeLength();
                for (std::size_t k = 0; k != n; ++k) {
                    real_t z = (k + 0.5) * edgeLength();
                    position_t r{x, y, z};
                    indices_t indices = std::make_tuple(i,j,k);
                    // Create cell without any particles.
                    cell_t cell{r, indices};                   
                    cells[i][j][k] = cell;
                }
            }
        }
        return std::make_tuple(cells, edgeLength, n);
    }
    
    template <typename P>
    static void 
    placeInCells_(boost::multi_array<Cell<P>, 3>& cells,
                  const length_t& edgeLength,
                  const bc_ptr_t& bc,
                  const std::vector<std::shared_ptr<P>>& particles,
                  const std::vector<std::shared_ptr<P>>& free,
                  const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {    
        /*
        for (auto p : particles) {
            auto r = p->position();
            r = bc->placeInside(r);
            std::size_t n0 = r[0] / edgeLength();
            std::size_t n1 = r[1] / edgeLength();
            std::size_t n2 = r[2] / edgeLength();
            auto& cell = cells[n0][n1][n2];
            auto& all = cell.all_;
            all.push_back(p);
        }
        */
        for (auto p : free) {
            auto r = p->position();
            r = bc->placeInside(r);
            std::size_t n0 = r[0] / edgeLength();
            std::size_t n1 = r[1] / edgeLength();
            std::size_t n2 = r[2] / edgeLength();
            auto& cell = cells[n0][n1][n2];
            auto& free = cell.free_;
            free.push_back(p);
        }
        for (auto g: groups) {
            auto r = g->position();
            r = bc->placeInside(r);
            std::size_t n0 = r[0] / edgeLength();
            std::size_t n1 = r[1] / edgeLength();
            std::size_t n2 = r[2] / edgeLength();
            auto& cell = cells[n0][n1][n2];
            auto& groups = cell.groups_;
            groups.push_back(g);
        }
    }
    
    template <typename P>
    static void 
    clearCells_(boost::multi_array<Cell<P>, 3>& cells,
                std::size_t size)
    {
        using cell_t = Cell<P>;
        
        for (auto iter = cells.origin(); 
                  iter != cells.origin() + cells.num_elements(); 
                ++iter) {
            cell_t& cell = *iter;
            cell.all_.clear();
            cell.free_.clear();
            cell.groups_.clear();
        }
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
            for (auto p : cell.free_) {
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
            for (auto pg: cell.groups_) {
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
            for (auto gj: cell.groups_) {
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
        using p_pair_list_t = typename CellLists<P>::p_pair_list_t;
        
        static boost::multi_array<Cell<P>, 3> cells{};
        static std::size_t size{0};
        static length_t edgeLength{0.0};
        static real_t rc2;
        
        // Set up if still required.
        static bool setup = false;
        if ( !setup ) {
            auto rc = rc_(box);
            auto created = makeCells_<P>(box, rc);
            cells = std::get<0>(created);
            edgeLength = std::get<1>(created);
            size = std::get<2>(created);
            rc2 = rc2_(box);
            setup = true;
        }
        
        // Determine full particle pair list.
        p_pair_list_t full{};
        clearCells_(cells, size);
        placeInCells_(cells, edgeLength, bc, all, free, groups);                
        for (auto iter = cells.origin(); 
                  iter != cells.origin() + cells.num_elements(); 
                ++iter) {
            const cell_t& cell = *iter;
            const auto& free = cell.free_;
            //const std::vector<typename CellLists<P>::p_ptr_t>& free = cell.free_;
            const auto& groups = cell.groups_;
            //const std::vector<typename CellLists<P>::pg_ptr_t>& groups = cell.groups_;
            auto neighbors = neighbors_(cell, size);
            for (std::size_t i = std::get<0>(neighbors); i != std::get<1>(neighbors); ++i) {
                for (std::size_t j = std::get<2>(neighbors); j != std::get<3>(neighbors); ++j) {
                    for (std::size_t k = std::get<4>(neighbors); k != std::get<5>(neighbors); ++k) {
                        cell_t& c = cells[i][j][k];
                        p_pair_list_t pl = forFree_<P>(free, c, bc, rc2);
                        full.insert(full.end(), pl.begin(), pl.end());
                        pl = forGroups_<P>(groups, c, bc, rc2);
                        full.insert(full.end(), pl.begin(), pl.end());
                    }
                }
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