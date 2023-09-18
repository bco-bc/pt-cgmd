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
     * Finds non-bonded particle pairs on the basis of distances between particles.
     */
    class DistancePairListGenerator : public pair_list_generator {
    public:

        /**
         * Constructor.
         * @param param Simulation parameters.
         * @param bc Boundary conditions.
         */
        DistancePairListGenerator(param_ptr_t param, bc_ptr_t bc);

        std::vector<PairList::p_pair_t>
        generate(const p_system_ptr_t& particleSystem) const override;

    private:

        param_ptr_t param_;
        bc_ptr_t bc_;
    };


}

#endif //SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP
