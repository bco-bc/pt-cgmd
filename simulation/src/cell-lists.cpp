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
#include <tuple>
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
                
    // Cell in grid. Holds free particles and groups.
    template <typename P>
    struct Cell {
        
        // Types.
        using p_ptr_t = typename CellLists<P>::p_ptr_t;
        using pg_ptr_t = typename CellLists<P>::pg_ptr_t;
        using indices_t = std::tuple<std::size_t, std::size_t, std::size_t>;
        
        // Constructors
        Cell() : r_{}, free_{}, groups_{}, indices_{} {}        
        Cell(const position_t& r, const indices_t& indices) : 
            r_{r}, free_{}, groups_{}, indices_{indices} {}
            
        position_t r_;                  // Cell position, center of cell.
        std::vector<p_ptr_t> free_;     // Free particles located in this cell.
        std::vector<pg_ptr_t> groups_;  // Particle groups located in this cell.
        indices_t indices_;             // Cell indices.
    };
    
    // Grid. Made off cells.
    template <typename P>
    struct Grid {

        // Constructors. Note: n is size per dimension.
        Grid() : cells_{} {}        
        Grid(std::size_t n) : 
            n_{n}, n3_{n * n * n}, cells_{n3_, Cell<P>{}}, edgeLength_{0.0} {}
        
        /**
         * Returns empty, but with position and indices of cells assigned.
         * @param box Box.
         * @param rc Cutoff distance.
         * @return Grid.
         */
        static Grid<P> make(const box_ptr_t& box,
                       const length_t& rc);
        
        // Getters.
        Cell<P>& operator () (std::size_t i, std::size_t j, std::size_t k);

        // Setters. Use cell indices.
        void set(const Cell<P>& cell);                
               
        // Returns location in cells_.
        std::size_t index(std::size_t i, std::size_t j, std::size_t k) const {
            return 9 * i + 3 * j + k;
        }
        
        /**
         * Returns neighboring cells.
         */
        std::vector<Cell<P>> neighbors(const Cell<P>& cell);
        
        /**
         * Clear free particles and groups from cell.s
         */
        void clear();
        
        /**
         * Place free particles and groups in cells.
         */
        void place(const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& particles,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups);
        
        std::size_t n_;
        std::size_t n3_;
        std::vector<Cell<P>> cells_;
        length_t edgeLength_;
    };
    
    template <typename P>
    Cell<P>& 
    Grid<P>::operator () (std::size_t i, std::size_t j, std::size_t k)
    {
        auto index = this->index(i, j, k);
        return cells_[index];
    }
    
    template <typename P>
    void 
    Grid<P>::set(const Cell<P>& cell)
    {
        auto i = std::get<0>(cell.indices_);
        auto j = std::get<1>(cell.indices_);
        auto k = std::get<2>(cell.indices_);
        auto index = this->index(i, j, k);
        assert(index < n3_);
        cells_[index] = cell;
    }
    
    template <typename P>
    std::vector<Cell<P>> 
    Grid<P>::neighbors(const Cell<P>& cell)
    {
        
        std::vector<Cell<P>> neighbors;

        const auto& indices = cell.indices_;       
        auto l = std::get<0>(indices);
        auto i_i = l > 0 ? l - 1 : l;
        auto i_f = l < n_ - 1 ? n_ : l;
        l = std::get<1>(indices);
        auto j_i = l > 0 ? l - 1 : l;
        auto j_f = l < n_ - 1 ? n_ : l; 
        l = std::get<2>(indices);
        auto k_i = l > 0 ? l - 1 : l;
        auto k_f = l < n_ - 1 ? n_ : l; 
        
        for ( std::size_t i = i_i; i != i_f; ++i) {
            for ( std::size_t j = j_i; j != j_f; ++j) {
                for ( std::size_t k = k_i; k != k_f; ++k) {
                    auto cell = this->operator ()(i, j, k);
                    neighbors.push_back(cell);
                }
            }
        }
        return std::move(neighbors);
    }
    
    /**
     * Makes the grid.
     * @param box Simulation box.
     * @param rc Cutoff distance.
     * @return Cells, cell edge length, and size per dimension.
     */
    template <typename P>
    Grid<P> 
    Grid<P>::make(const box_ptr_t& box, 
                  const length_t& rc)
    {
        using cell_t = Cell<P>;
        using indices_t = typename Cell<P>::indices_t;
        
        auto edgeLength = 0.5 * rc;
        std::size_t n = box->edgeLength() / edgeLength();
        Grid<P> grid{n};
        grid.edgeLength_ = edgeLength;
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
                    grid.set(cell);
                }
            }
        }
        return grid;
    }
    
    template <typename P>
    void 
    Grid<P>::clear()
    {
        for (auto iter = cells_.begin(); iter != cells_.end(); ++iter) {
            auto& cell = *iter;
            cell.free_.clear();
            cell.groups_.clear();
        }
    }
    
    
    template <typename P>
    void 
    Grid<P>::place(const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& particles,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {    
        for (auto p : free) {
            auto r = p->position();
            r = bc->placeInside(r);
            std::size_t i = r[0] / edgeLength_();
            std::size_t j = r[1] / edgeLength_();
            std::size_t k = r[2] / edgeLength_();
            auto& cell = this->operator()(i, j, k);
            auto& free = cell.free_;
            free.push_back(p);
        }
        for (auto g: groups) {
            auto r = g->position();
            r = bc->placeInside(r);
            std::size_t i = r[0] / edgeLength_();
            std::size_t j = r[1] / edgeLength_();
            std::size_t k = r[2] / edgeLength_();
            auto& cell = this->operator()(i, j, k);
            auto& groups = cell.groups_;
            groups.push_back(g);
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
        grid.place(bc, all, free, groups);
        auto& cells = grid.cells_;
        for (auto iter = cells.begin(); iter != cells.end(); ++iter) {
            const cell_t& cell = *iter;
            const auto& free = cell.free_;
            const auto& groups = cell.groups_;
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