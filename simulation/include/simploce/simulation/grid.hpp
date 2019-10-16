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
     * cell lists. The grid and its elements are always cubically shaped. 
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
         * Type of cell location in grid.
         */
        using location_t = typename cell_t::location_t;

        /**
         * Constructor. Empty grid, no cells.
         */
        Grid();
        
        /**
         * Creates new grid, with location and position assigned to individual 
         * cells. The number of cells in one dimension 
         * is determined from the ratio of the box size and element side length.
         * @param box Box.
         * @param sideLength Requested element side length. The actual side length 
         * may be slightly different.
         * @return Grid.
         */
        static Grid<P> make(const box_ptr_t& box,
                            const length_t& sideLength);
        
        /**
         * Returns total number of cells.
         * @return Number of cells.
         */
        std::size_t numberOfCells() const { return ncells_; }
        
        /**
         * Returns an individual cell.
         * @param location Location in grid.
         * @return Cell, is modifiable.
         */
        cell_t& operator () (const location_t& location);

        /**
         * Replaces an existing cell at the same location in the grid as 
         * the given cell.
         * @param cell Cell.
         */
        void replace(const Cell<P>& cell);                
               
        /**
         * Returns neighboring cells around a cell with given location.
         * @param location Location in grid.
         */
        std::vector<cell_t> neighbors(const location_t& location);
        
        /**
         * Clears free particles and groups from cells.
         */
        void clear();
        
        /**
         * Places free particles and particles groups in cells.
         * @param bc Boundary condition.
         * @param free Free particles.
         * @param groups Particle groups.
         */
        void place(const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups);
        
        /**
         * Returns all cells.
         * @return Cells.
         */
        const std::vector<cell_t>& cells() { return cells_; }
        
    private:
        
        /**
         * Constructor.
         * @param n Number of grid elements in one dimension, and total number 
         * of grid elements is therefore n * n * n. The three element indices
         * jointly specifying location of the cell in the grid, should be in 
         * the set {0, 1, n-1}.
         * @param sideLength. Used to determine the position of cells in the grid.
         */
        Grid(std::size_t n, const length_t& sideLength);
        
        std::size_t index_(const location_t& location) const;
        
        static location_t location_(std::size_t i, std::size_t j, std::size_t k) {
            return std::make_tuple(i, j, k);
        }
        
        std::size_t ncells_1_;  // Number of cells in one direction.
        std::size_t ncells_;    
        std::vector<cell_t> cells_;
        length_t sideLength_;
        
    };
    
    template <typename P>
    Grid<P>::Grid() :
        ncells_1_{0}, ncells_{0}, cells_{}, sideLength_{0.0}
    {        
    }
    
    template <typename P>
    Grid<P>::Grid(std::size_t n, const length_t& sideLength) : 
        ncells_1_{n}, ncells_{n * n * n}, cells_{}, sideLength_{sideLength} 
    {
        for ( std::size_t i = 0; i != n; ++i) {
            real_t x = (i + 0.5) * sideLength_();
            for ( std::size_t j = 0; j != n; ++j) {
                real_t y = (j + 0.5) * sideLength_();
                for ( std::size_t k = 0; k != n; ++k) {
                    real_t z = (k + 0.5) * sideLength_();
                    position_t r{x, y, z};
                    location_t location = this->location_(i, j, k);
                    cell_t cell{r, location};
                    cells_.push_back(cell);
                }
            }
        }
    }
        
    template <typename P>
    typename Grid<P>::cell_t& 
    Grid<P>::operator () (const location_t& location)
    {
        auto index = this->index_(location);
        return cells_[index];
    }
    
    template <typename P>
    void 
    Grid<P>::replace(const cell_t& cell)
    {
        auto index = this->index_(cell.location());
        cells_[index] = cell;
    }
    
    template <typename P>
    std::vector<typename Grid<P>::cell_t> 
    Grid<P>::neighbors(const location_t& location)
    {        
        std::vector<cell_t> neighbors{};        
        std::size_t n = ncells_1_;

        // Specify lower and upper limits.
        auto l = std::get<0>(location);
        auto i_i = l > 0 ? l - 1 : l;
        auto i_f = l < n - 1 ? n : l;
        l = std::get<1>(location);
        auto j_i = l > 0 ? l - 1 : l;
        auto j_f = l < n - 1 ? n : l; 
        l = std::get<2>(location);
        auto k_i = l > 0 ? l - 1 : l;
        auto k_f = l < n - 1 ? n : l; 
        
        // Find associated cells.
        for ( std::size_t i = i_i; i != i_f; ++i) {
            for ( std::size_t j = j_i; j != j_f; ++j) {
                for ( std::size_t k = k_i; k != k_f; ++k) {
                    auto location = this->location_(i, j, k);
                    const auto& cell = this->operator ()(location);
                    neighbors.push_back(cell);
                }
            }
        }
        return std::move(neighbors);
    }
    
    template <typename P>
    Grid<P> 
    Grid<P>::make(const box_ptr_t& box, 
                  const length_t& sideLength)
    {
        std::size_t n = util::nint(box->edgeLength() / sideLength());
        real_t nreal = n;
        real_t edgeLength = box->edgeLength() / nreal;
        return std::move(Grid<P>{n, edgeLength});
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
        auto el = sideLength_();
        for (auto p : free) {
            auto r = p->position();
            r = bc->placeInside(r);
            std::size_t i = r[0] / el;
            std::size_t j = r[1] / el;
            std::size_t k = r[2] / el;
            auto location = this->location_(i, j, k);
            auto& cell = this->operator()(location);
            auto& free = cell.free_;
            free.push_back(p);
        }
        for (auto g : groups) {
            auto r = g->position();
            r = bc->placeInside(r);
            std::size_t i = r[0] / sideLength_();
            std::size_t j = r[1] / sideLength_();
            std::size_t k = r[2] / sideLength_();
            auto location = this->location_(i, j, k);
            std::cout << "i, j, k, index, ncells_ "
                      << i << ' ' << j << ' ' << k
                      << ' ' << this->index_(location) << std::endl;
            auto& cell = this->operator()(location);
            cell.groups_.push_back(g);
        }
    }
    
    template <typename P>
    std::size_t 
    Grid<P>::index_(const location_t& location) const {
        auto i = std::get<0>(location);
        auto j = std::get<1>(location);
        auto k = std::get<2>(location);
        std::size_t index = 9 * i + 3 * j + k;
        assert(index < ncells_);
        return index;
    }
    
}

#endif /* GRID_HPP */

