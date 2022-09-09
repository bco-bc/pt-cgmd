/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 15, 2019, 3:26 PM
 */

#ifndef GRID_HPP
#define GRID_HPP

#include "simploce/types/s-types.hpp"
#include <map>

namespace simploce {
    
    /**
     * Grid, consists of cells.
     */
    class Grid {
    public:

        /**
         * Constructor
         * @param box Simulation box.
         * @param bc Boundary conditions.
         * @param cutoff Cutoff distance for cell pairs.
         */
        Grid(box_ptr_t box, bc_ptr_t bc, dist_t cutoff);
        
         /**
         * Places free particles and particles groups in cells.
         * @param bc Boundary condition.
         * @param particles Particles involved in non-bonded interactions.
         */
        void 
        assignParticlesToCells(const bc_ptr_t& bc,
                               const std::vector<p_ptr_t>& particles);

        /**
         * Return neighbors of each cell.
         * @return Neighbors. The key of the std::map map corresponds to the location of any given cell
         * in the grid.
         */
        const std::map<std::string, std::vector<cell_ptr_t>>& neighbors() const;

        /**
         * Returns a cell.
         * @param location Location of the cell.
         * @return Cell.
         */
        cell_ptr_t findCell(const std::string& location) const;

        /**
         * Returns cells
         * @return Cells.
         */
        const std::vector<cell_ptr_t>& cells() const;

        /**
         * Returns cutoff distance for particle pairs.
         * @return Cutoff distance.
         */
        dist_t cutoff() const;

    private:

        // Clears assigned particles from cells.
        void clear_();

        box_ptr_t box_;
        bc_ptr_t bc_;
        dist_t cutoff_;

        std::map<std::string, cell_ptr_t> cells_{};
        int nx_{0};
        real_t spacingX_{0.0};
        int ny_{0};
        real_t spacingY_{0.0};
        int nz_{0};
        real_t spacingZ_{0.0};
        std::map<std::string, std::vector<cell_ptr_t>> neighbors_{};
    };
    

}

#endif /* GRID_HPP */

