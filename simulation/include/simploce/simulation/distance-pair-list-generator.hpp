/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#ifndef SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP
#define SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP

#include "pair-list-generator.hpp"

namespace simploce {

    /**
     * Generates pair lists on the basis of distances between particle and particle groups.e.
     */
    class DistancePairListGenerator : public pair_lists_generator {
    public:

        /**
         * Constructor.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        DistancePairListGenerator(param_ptr_t param, bc_ptr_t bc);

        PairLists generate(const p_system_ptr_t& particleSystem) const override;

    private:

        param_ptr_t param_;
        box_ptr_t box_;
        bc_ptr_t bc_;
    };


}

#endif //SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP
