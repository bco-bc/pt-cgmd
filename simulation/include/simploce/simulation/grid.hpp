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
         * Grid pointer type.
         */
        using grid_ptr_t = std::shared_ptr<Grid<P>>;
    
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<ParticleGroup<P>>;
        
        /**
         * Cell type.
         */
        using cell_t = Cell<P>;
        
        /**
         * Type of cell location in grid.
         */
        using location_t = typename cell_t::location_t;

        /**
         * Creates new grid, with location and position assigned to individual 
         * cells. The number of cells in one dimension is determined from the 
         * ratio of the box size and element side length.
         * @param box Box.
         * @param sideLength Requested element side length. The actual side length 
         * of the new grid may be different.
         * @return Grid.
         */
        static grid_ptr_t 
        make(const box_ptr_t& box,
             const length_t& sideLength);
        
        /**
         * Returns total number of cells.
         * @return Number of cells.
         */
        std::size_t 
        numberOfCells() const { return ncells_; }
        
        /**
         * Returns the number of cell along one axis or dimension.
         * @return Number of cells.
         */
        std::size_t
        numberOfCellPerDimension() const { return ncells1D_; }
        
        /**
         * Returns an individual cell.
         * @param location Location in grid.
         * @return Cell, is modifiable.
         */
        cell_t& 
        operator () (const location_t& location);

        /**
         * Replaces an existing cell at the same location in the grid as 
         * the given cell.
         * @param cell Cell.
         */
        void 
        replace(const cell_t& cell);                
               
        /**
         * Returns neighboring cells around a cell with given location.
         * @param location Location in grid.
         */
        std::vector<cell_t> 
        neighbors(const location_t& location);
        
        /**
         * Clears cells from free particles and particle groups.
         */
        void 
        clear();
        
        /**
         * Places free particles and particles groups in cells.
         * @param bc Boundary condition.
         * @param free Free particles.
         * @param groups Particle groups.
         */
        void 
        place(const bc_ptr_t& bc,
              const std::vector<p_ptr_t>& free,
              const std::vector<pg_ptr_t>& groups);
        
        /**
         * Returns all cells.
         * @return Cells.
         */
        const std::vector<cell_t>& cells() { return cells_; }
        
        /**
         * Returns the side length of cells in this grid.
         * @return Side length.
         */
        length_t sideLength() const { return sideLength_; }
        
    private:
                
        Grid(std::size_t n, 
            const length_t& sideLength);
        
        std::size_t 
        index_(const location_t& location) const;
        
        static location_t location_(std::size_t i, std::size_t j, std::size_t k) {
            return std::make_tuple(i, j, k);
        }
        
        std::size_t ncells1D_;       // Number of cells in one direction.
        std::size_t ncells_;    
        std::vector<cell_t> cells_;
        length_t sideLength_;
        
    };
    
    template <typename P>
    Grid<P>::Grid(std::size_t n, 
                  const length_t& sideLength) : 
        ncells1D_{n}, ncells_{n * n * n}, cells_{}, sideLength_{sideLength} 
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
        //std::clog << "index: " << this->index_(location) << std::endl;
        std::vector<cell_t> neighbors{};        
        int n = ncells1D_;
        
        // Specify lower and upper limits in each direction.
        std::size_t lx = std::get<0>(location);  // 0 <= lx < n
        int start = lx - 1;
        std::size_t iStart = (start < 0 ? 0 : start);
        int end = lx + 2;
        std::size_t iEnd = (end > n ? n : end);
        
        std::size_t ly = std::get<1>(location);  // 0 <= ly < n
        start = ly - 1;
        std::size_t jStart = (start < 0 ? 0 : start);
        end = ly + 2;
        std::size_t jEnd = (end > n ? n : end);
        
        std::size_t lz = std::get<2>(location);  // 0 <= lz < n
        start = lz - 1;
        std::size_t kStart = (start < 0 ? 0 : start);
        end = lz + 2;
        std::size_t kEnd = (end > n ? n : end);
        
        // Select neighboring cells.
        for ( std::size_t i = iStart; i != iEnd; ++i) {
            for ( std::size_t j = jStart; j != jEnd; ++j) {
                for ( std::size_t k = kStart; k != kEnd; ++k) {
                    // Exclude itself
                    if ( !(i == lx && j == ly && k == lz) ) {
                        auto location = this->location_(i, j, k);
                        const auto& cell = this->operator ()(location);
                        neighbors.push_back(cell);
                    }
                }
            }
        }
        return std::move(neighbors);
    }
    
    template <typename P>
    typename Grid<P>::grid_ptr_t
    Grid<P>::make(const box_ptr_t& box, 
                  const length_t& sideLength)
    {
        // Adapts side lengths.
        std::size_t n = box->edgeLength() / sideLength();
        real_t edgeLength = box->edgeLength() / n;
        return grid_ptr_t{(new Grid<P>(n, edgeLength))};
    }
    
    template <typename P>
    void 
    Grid<P>::clear()
    {
        for (auto& cell : cells_) {
            cell.clear();
        }
    }
    
    
    template <typename P>
    void 
    Grid<P>::place(const bc_ptr_t& bc,
                   const std::vector<std::shared_ptr<P>>& free,
                   const std::vector<std::shared_ptr<ParticleGroup<P>>>& groups)
    {    
        for (auto p : free) {
            auto r = bc->placeInside(p->position());
            std::size_t i = r[0] / sideLength_();
            std::size_t j = r[1] / sideLength_();
            std::size_t k = r[2] / sideLength_();
            auto location = this->location_(i, j, k);
            auto& cell = this->operator()(location);
            auto& free = cell.free_;
            free.push_back(p);
        }
        for (auto g : groups) {
            auto r = bc->placeInside(g->position());
            std::size_t i = r[0] / sideLength_();
            std::size_t j = r[1] / sideLength_();
            std::size_t k = r[2] / sideLength_();
            auto location = this->location_(i, j, k);                        
            auto& cell = this->operator()(location);
            cell.groups_.push_back(g);
        }
    }
    
    template <typename P>
    std::size_t 
    Grid<P>::index_(const location_t& location) const 
    {
        auto i = std::get<0>(location);
        auto j = std::get<1>(location);
        auto k = std::get<2>(location);
        std::size_t index = ncells1D_ * ncells1D_ * i + ncells1D_ * j + k;
        assert(index < ncells_);
        return index;
    }
    
}

#endif /* GRID_HPP */

