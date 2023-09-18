/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11 October 2019, 16:29
 */

#ifndef CELL_PAIR_LIST_GENERATOR_HPP
#define CELL_PAIR_LIST_GENERATOR_HPP

#include "pair-list-generator.hpp"
#include "simploce/types/s-types.hpp"

namespace simploce {
    
    /**
     * Finds non-bonded particle pairs based on the location of particles in cells.
     * @see <a href="https://en.wikipedia.org/wiki/Cell_lists">Cell lists at Wikipedia</a>
     */
    class CellPairListGenerator : public pair_list_generator {
    public:

        /**
         * Constructor.
         * @param param Simulation parameters.
         * @param bc Boundary conditions.
         */
        CellPairListGenerator(param_ptr_t param, bc_ptr_t bc);
        
        std::vector<PairList::p_pair_t>
        generate(const p_system_ptr_t& particleSystem) const override;
        
    private:
                        
        param_ptr_t param_;
        bc_ptr_t bc_;
    };

}

#endif /* CELL_PAIR_LIST_GENERATOR_HPP */

