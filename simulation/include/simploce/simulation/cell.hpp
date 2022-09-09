/*
 * File:   findCell.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 15, 2019, 3:10 PM
 */

#ifndef CELL_HPP
#define CELL_HPP

#include "simploce/types/s-types.hpp"
#include "simploce/particle/particle-group.hpp"
#include <vector>
#include <tuple>
#include <memory>
#include <string>

namespace simploce {
    
    /**
     * An element of a grid. Holds a list of free particles and particle groups
     * located in this findCell.
     */
    class Cell {
    public:
        
        /*
         * Indices specifying location in grid as multiples of the grid spacing (length cell edge) in 3
         * directions.
         */
        using location_t = std::tuple<std::size_t, std::size_t, std::size_t>;

        /**
         * Returns location as std::string.
         * @param location
         * @return std::string. Of the form 'i-j-k'.
         */
        static std::string stringLocation(const location_t& location);
        
        /**
         * Constructor. Creates cell with location in grid. No particles assigned.
         * @param location Location in grid.
         * @param r Position in cartesian coordinate system.
         */
        Cell(location_t location, position_t r);

        // Noncopyable
        Cell(const Cell&) = delete;
        Cell& operator = (const Cell&) = delete;

        /**
         * Equality operator.
         * @param cell Cell.
         * @return Result.
         */
        bool operator == (const Cell& cell) const;

        /**
         * Assigns given particle to this cell.
         * @param particle Particle.
         */
        void assign(const p_ptr_t& particle);

        /**
         * Returns all particles currently located in/assigned to this cell.
         * @return Particles. Unmodifiable.
         */
        const std::vector<p_ptr_t>& particles() const;
        
        /**
         * Returns location in a grid.
         * @return Location.
         */
        location_t location() const;

        /**
         * Returns position of cell.
         * @return Position (cartesian coordinates)
         */
        position_t position() const;

        /**
         * Returns location as a string.
         * @return String, of the form "i-j-k".
         */
        std::string stringLocation() const;
        
    private:

        friend class Grid;

        /**
         * Clears all particles assigned to this cell.
         */
        void clear_();

        std::vector<p_ptr_t> particles_;
        location_t location_;
        position_t r_;
        std::string stringLocation_;
    };


    
}

#endif /* CELL_HPP */

