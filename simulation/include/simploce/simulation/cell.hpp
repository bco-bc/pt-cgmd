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
 * File:   cell.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 15, 2019, 3:10 PM
 */

#ifndef CELL_HPP
#define CELL_HPP

#include "stypes.hpp"
#include "simploce/particle/particle-group.hpp"
#include <vector>
#include <tuple>
#include <memory>

namespace simploce {
    
    /**
     * An element of a grid. Holds free particles and groups.
     * @param P Particle type.
     */
    template <typename P>
    class Cell {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = std::shared_ptr<ParticleGroup<P>>;
        
        /**
         * Indices. Specifies location in grid.
         */
        using indices_t = std::tuple<std::size_t, std::size_t, std::size_t>;
        
        /**
         * Constructor. Empty cell.
         */
        Cell() : r_{}, free_{}, groups_{}, indices_{} {}        
        
        /**
         * Constructor. Creates cell with position and location in grid.
         * @param r
         * @param indices
         */
        Cell(const position_t& r, const indices_t& indices) : 
            r_{r}, free_{}, groups_{}, indices_{indices} {}
            
        /**
         * Clears this cell from free particles and particle groups.
         */
        void clear() { free_.clear(); groups_.clear(); }
        
        /**
         * Returns free particles.
         * @return Free particles.
         */
        const std::vector<p_ptr_t>& free() const { return free_; }
        
        /**
         * Returns particle groups.
         * @return Particle groups.
         */
        const std::vector<pg_ptr_t>& groups() const { return groups_; }
            
    private:
        
        template <typename PP>
        friend class Grid;
        
        position_t r_;                  // Cell position, center of cell.
        std::vector<p_ptr_t> free_;     // Free particles located in this cell.
        std::vector<pg_ptr_t> groups_;  // Particle groups located in this cell.
        indices_t indices_;             // Cell indices.
    };
    
}

#endif /* CELL_HPP */

