/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/8/21.
 */

#ifndef SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP
#define SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP

#include "pair-list-generator.hpp"
#include "pair-lists.hpp"
#include "s-types.hpp"

namespace simploce {


    /**
     * Creates atom pair lists.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @param all All atoms.
     * @param free Free atoms.
     * @param groups Atom groups.
     * @return Atom pair lists.
     */
    PairLists<Atom>
    createPairLists(const box_ptr_t& box,
                    const bc_ptr_t& bc,
                    const std::vector<atom_ptr_t>& all,
                    const std::vector<atom_ptr_t>& free,
                    const std::vector<atom_group_ptr_t>& groups);

    /**
     * Creates bead pair lists.
     * @param box Simulation box.
     * @param bc Boundary condition.
     * @param all All beads.
     * @param free Free beads.
     * @param groups Bead groups.
     * @return Bead pair lists.
     */
    PairLists<Bead>
    createPairLists(const box_ptr_t& box,
                    const bc_ptr_t& bc,
                    const std::vector<bead_ptr_t>& all,
                    const std::vector<bead_ptr_t>& free,
                    const std::vector<bead_group_ptr_t>& groups);

    /**
     * Generates pair lists on the basis of distances between particle and particle groups.
     * @tparam P Particle type.
     */
    template <typename P>
    class DistancePairListGenerator : public pair_lists_generator<P> {
    public:

        using p_ptr_t = typename pair_lists_generator<P>::p_ptr_t;

        using pg_ptr_t = typename pair_lists_generator<P>::pg_ptr_t;

        DistancePairListGenerator(const box_ptr_t &box,
                                  const bc_ptr_t &bc);

        ~DistancePairListGenerator();

        PairLists<P>
        generate(const std::vector<p_ptr_t>& all,
                 const std::vector<p_ptr_t>& free,
                 const std::vector<pg_ptr_t>& groups) const override;
    private:

        box_ptr_t box_;
        bc_ptr_t bc_;
    };

    template <typename P>
    DistancePairListGenerator<P>::DistancePairListGenerator(const box_ptr_t &box,
                                                            const bc_ptr_t &bc) :
        pair_lists_generator<P>{}, box_{box}, bc_{bc} {
    }

    template <typename P>
    DistancePairListGenerator<P>::~DistancePairListGenerator() = default;

    template <typename P>
    PairLists<P>
    DistancePairListGenerator<P>::generate(const std::vector<p_ptr_t>& all,
                                           const std::vector<p_ptr_t>& free,
                                           const std::vector<pg_ptr_t>& groups) const {
         return createPairLists(box_, bc_, all, free, groups);
    }

}

#endif //SIMULATION_DISTANCE_PAIR_LIST_GENERATOR_HPP
