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
     * Creates atom pair lists.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @param all All atoms.
     * @param free Free atoms.
     * @param groups Atom groups.
     * @return Atom pair lists.
     *
    PairLists
    createPairLists(const box_ptr_t& box,
                    const bc_ptr_t& bc,
                    const std::vector<p_ptr_t>& all,
                    const std::vector<p_ptr_t>& free,
                    const std::vector<pg_ptr_t>& groups);
                    */

    /**
     * Generates pair lists on the basis of distances between particle and particle groups.
     * @tparam P Particle type.
     */
    class DistancePairListGenerator : public pair_lists_generator {
    public:

        DistancePairListGenerator(box_ptr_t box,
                                  bc_ptr_t bc);

        ~DistancePairListGenerator();

        PairLists
        generate(const std::vector<p_ptr_t>& all,
                 const std::vector<p_ptr_t>& free,
                 const std::vector<pg_ptr_t>& groups) const override;
    private:

        box_ptr_t box_;
        bc_ptr_t bc_;
    };


}

#endif //SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP
