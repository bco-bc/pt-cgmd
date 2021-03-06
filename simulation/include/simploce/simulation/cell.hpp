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
        using location_t = std::tuple<std::size_t, std::size_t, std::size_t>;
        
        /**
         * Constructor. Creates cell with position and location in grid. No
         * free particles and particles groups.
         * @param r Position in Cartesian coordinates system.
         * @param indices Location in grid.
         */
        Cell(const position_t& r, const location_t& location) : 
            r_{r}, free_{}, groups_{}, location_{location} {} 
            
        /**
         * Clears this cell from free particles and particle groups.
         */
        void clear() { free_.clear(); groups_.clear(); }
        
        /**
         * Returns position in Cartesian coordinate system.
         * @return Position.
         */
        position_t position() const { return r_; }
        
        /**
         * Returns free particles currently located in this cell.
         * @return Free particles.
         */
        const std::vector<p_ptr_t>& free() const { return free_; }
        
        /**
         * Returns particle groups currently located in this cell.
         * @return Particle groups.
         */
        const std::vector<pg_ptr_t>& groups() const { return groups_; }
        
        /**
         * Returns cell location in grid.
         * @return Location.
         */
        location_t location() const { return location_; }
            
    private:
        
        template <typename PP>
        friend class Grid;
        
        position_t r_;
        std::vector<p_ptr_t> free_;
        std::vector<pg_ptr_t> groups_;
        location_t location_;
    };
    
}

#endif /* CELL_HPP */

