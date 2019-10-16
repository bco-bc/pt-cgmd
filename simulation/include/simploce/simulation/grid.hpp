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
 * File:   grid.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 15, 2019, 3:26 PM
 */

#ifndef GRID_HPP
#define GRID_HPP

#include "stypes.hpp"
#include "cell.hpp"
#include "bc.hpp"
#include "simploce/util/util.hpp"
#include <vector>
#include <memory>
#include <utility>

namespace simploce {
    
    /**
     * Grid, consists of cells. For creating particle based pair lists based on 
     * cell lists.
     * @param P Particle type.
     */
    template <typename P>
    class Grid {
    public:
    
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Cell type.
         */
        using cell_t = Cell<P>;

        /**
         * Constructor. Empty grid.
         */
        Grid() : cells_{} {}
        
        /**
         * Constructor.
         * @param n Size of grid in one direction. Total number of grid elements
         * is n * n * n.
         */
        Grid(std::size_t n);
        
        /**
         * Returns empty, but with position and cell indices assigned to individual 
         * cells.
         * @param box Box.
         * @param rcutoff Cutoff distance.
         * @return Grid.
         */
        static Grid<P> make(const box_ptr_t& box,
                            const length_t& rcutoff);
        
        /**
         * Returns cell.
         * @param i Index in x direction.
         * @param j Index in y direction.
         * @param k Index in z direction.
         * @return Cell, is modifiable.
         */
        cell_t& operator () (std::size_t i, 
                             std::size_t j, 
                             std::size_t k);

        /**
         * Replaces cell.
         * @param cell Cell that replaces an existing cell.
         */
        void replace(const Cell<P>& cell);                
               
        /** 
         * Returns location in grid..
         * @param i Index in x-direction.
         * @param i Index in x-direction. 
         * @param i Index in x-direction.
         * @return Index.
         */
        std::size_t index(std::size_t i, 
                          std::size_t j, 
                          std::size_t k) const { return 9 * i + 3 * j + k; }
        
        /**
         * Returns neighboring cells around given cell.
         * @param cell Cell. Should be part of this grid.
         */
        std::vector<cell_t> neighbors(const cell_t& cell);
        
        /**
         * Clear free particles and groups from cells.
         */
        void clear();
        
        /**
         * Place free particles and groups in cells.
         * @param bc Boundary condition.
         * @param free Free particles.
         * @param groups Particle groups.
         */
        void place(const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups);
        
        /**
         * Returns cells.
         * @return Cells.
         */
        const std::vector<cell_t>& cells() { return cells_; }
        
    private:
        
        std::size_t ncells_;
        std::vector<cell_t> cells_;
        length_t edgeLength_;
        
    };
    
    template <typename P>
    Grid<P>::Grid(std::size_t n) : 
        ncells_{n * n * n}, cells_(ncells_, Cell<P>{}), edgeLength_{0.0} 
    {
    }
    
    template <typename P>
    typename Grid<P>::cell_t& 
    Grid<P>::operator () (std::size_t i, 
                          std::size_t j, 
                          std::size_t k)
    {
        auto index = this->index(i, j, k);
        return cells_[index];
    }
    
    template <typename P>
    void 
    Grid<P>::replace(const cell_t& cell)
    {
        auto i = std::get<0>(cell.indices_);
        auto j = std::get<1>(cell.indices_);
        auto k = std::get<2>(cell.indices_);
        auto index = this->index(i, j, k);
        assert(index < ncells_);
        cells_[index] = cell;
    }
    
    template <typename P>
    std::vector<typename Grid<P>::cell_t> 
    Grid<P>::neighbors(const cell_t& cell)
    {
        
        std::vector<cell_t> neighbors{};
        
        std::size_t n = std::pow(real_t(ncells_), 1.0/3.0);

        const auto& indices = cell.indices_;       
        auto l = std::get<0>(indices);
        auto i_i = l > 0 ? l - 1 : l;
        auto i_f = l < n - 1 ? n : l;
        l = std::get<1>(indices);
        auto j_i = l > 0 ? l - 1 : l;
        auto j_f = l < n - 1 ? n : l; 
        l = std::get<2>(indices);
        auto k_i = l > 0 ? l - 1 : l;
        auto k_f = l < n - 1 ? n : l; 
        
        for ( std::size_t i = i_i; i != i_f; ++i) {
            for ( std::size_t j = j_i; j != j_f; ++j) {
                for ( std::size_t k = k_i; k != k_f; ++k) {
                    auto& cell = this->operator ()(i, j, k);
                    neighbors.push_back(cell);
                }
            }
        }
        return std::move(neighbors);
    }
    
    /**
     * Returns new grid.
     * @param box Simulation box.
     * @param rc Cutoff distance for non-bonded interactions.
     * @return Grid.
     */
    template <typename P>
    Grid<P> 
    Grid<P>::make(const box_ptr_t& box, 
                  const length_t& rc)
    {
        using indices_t = typename Cell<P>::indices_t;
        
        std::size_t n = util::nint(box->edgeLength() / rc());
        real_t edgeLength = box->edgeLength() / n;
        Grid<P> grid{n};
        grid.edgeLength_ = edgeLength;
        for (std::size_t i = 0; i != n; ++i) {
            real_t x = (i + 0.5) * edgeLength;
            for (std::size_t j = 0; j != n; ++j) {
                real_t y = (j + 0.5) * edgeLength;
                for (std::size_t k = 0; k != n; ++k) {
                    real_t z = (k + 0.5) * edgeLength;
                    position_t r{x, y, z};
                    indices_t indices = std::make_tuple(i,j,k);
                    // Create cell without any particles.
                    cell_t cell{r, indices};                   
                    grid.replace(cell);
                }
            }
        }
        return std::move(grid);
    }
    
    template <typename P>
    void 
    Grid<P>::clear()
    {
        for (auto iter = cells_.begin(); iter != cells_.end(); ++iter) {
            auto& cell = *iter;
            cell.clear();
        }
    }
    
    
    template <typename P>
    void 
    Grid<P>::place(const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {    
        auto el = edgeLength_();
        for (auto p : free) {
            auto r = p->position();
            r = bc->placeInside(r);
            std::size_t i = r[0] / el;
            std::size_t j = r[1] / el;
            std::size_t k = r[2] / el;
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
            std::cout << "i,j,k,index,ncells_ "
                      << i << ' ' << j << ' ' << k
                      << ' ' << this->index(i, j, k) << std::endl;
            auto& cell = this->operator()(i, j, k);
            auto& groups = cell.groups_;
            groups.push_back(g);
        }
    }
    
    
}

#endif /* GRID_HPP */

