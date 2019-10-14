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
        for (std::size_t i = 0; i != size; ++i) {
            for (std::size_t j = 0; j != size; ++j) {
                for (std::size_t k = 0; k != size; ++k) {
                    auto& cell = cells[i][j][k];
                    cell.all_.clear();
                    cell.free_.clear();
                    cell.groups_.clear();
                }
            }
        }
    }
    
    template <typename P>
    static std::vector<typename CellLists<P>::p_pair_list_t>
    forFree(const boost::multi_array<Cell<P>, 3>& cells,
            std::size_t size,
            const bc_ptr_t& bc,
            real_t rc2)
    {
        using cell_t = Cell<P>;
        
        std::vector<typename CellLists<P>::p_pair_list_t> pairList{};
        
        for (auto iter = cells.begin(); iter != cells.end(); ++iter) {
            const cell_t& cell = *iter;
            auto neighbors = neighbors_(cell, size);
            for (std::size_t i = std::get<0>(neighbors); i != std::get<1>(neighbors); ++i) {
                for (std::size_t j = std::get<2>(neighbors); j != std::get<3>(neighbors); ++j) {
                    for (std::size_t k = std::get<4>(neighbors); k != std::get<5>(neighbors); ++k) {
                        
                    }
                }
            }
        }
        return pairList;
    }

    template <typename P>
    static std::vector<typename CellLists<P>::p_pair_list_t>
    forFreeAndGroups()
    {
        std::vector<typename CellLists<P>::p_pair_list_t> pairList{};
        
        return pairList;        
    }
    
    template <typename P>
    static std::vector<typename CellLists<P>::p_pair_list_t>
    forGroups()
    {
        std::vector<typename CellLists<P>::p_pair_list_t> pairList{};
        
        return pairList;        
    }
    
    template <typename P>
    static std::vector<typename CellLists<P>::p_pair_list_t>
    makePairLists_(const box_ptr_t& box,
                   const bc_ptr_t& bc,
                   const std::vector<typename CellLists<P>::p_ptr_t>& all,
                   const std::vector<typename CellLists<P>::p_ptr_t>& free,
                   const std::vector<typename CellLists<P>::pg_ptr_t>& groups)
    {
        using p_pair_list_t = typename CellLists<P>::p_pair_list_t;
        
        static boost::multi_array<Cell<P>, 3> cells{};
        static std::size_t size{0};
        static length_t edgeLength{0.0};
        
        static bool setup = false;
        if ( !setup ) {
            auto rc = rc_(box);
            auto created = makeCells_<P>(box, rc);
            cells = std::get<0>(created);
            edgeLength = std::get<1>(created);
            size = std::get<2>(created);
            setup = true;
        }
        
        clearCells_(cells, size);
        placeInCells_(cells, edgeLength, bc, all, free, groups);
        
        
        return std::vector<p_pair_list_t>{};
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
        return makePairLists_<Bead>(box_, bc_, all, free, groups);
    }
}