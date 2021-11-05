/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 2, 2019, 4:13 PM
 */

#ifndef CG_FORCEFIELD_HPP
#define CG_FORCEFIELD_HPP

#include "forcefield.hpp"
#include "pair-lists.hpp"
#include "s-types.hpp"
#include <utility>
#include <vector>

namespace simploce {
    
    /**
     * Interface for coarse grained force fields.
     */
    struct CoarseGrainedForceField : public ForceField {
                        
        virtual ~CoarseGrainedForceField() {}
        
        /**
         * Computes forces due to -all- interactions on beads. Updates/adds all forces 
         * acting on beads.
         * @param all All beads.
         * @param free Free beads.
         * @param groups All bead groups.
         * @param pairLists Pair lists.
         * @return Bonded and non-bonded potential energy.
         */
        virtual std::pair<energy_t, energy_t> 
        interact(const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups,
                 const PairLists<Bead>& pairLists) = 0;
        
        /**
         * Computes forces due to bonded interactions on beads. Updates/adds forces 
         * bonded acting on beads.
         * @param all All beads.
         * @param free Free beads.
         * @param groups  All bead groups.
         * @param pairLists Pair lists.
         * @return Potential energy for bonded interactions.
         */
        virtual energy_t 
        bonded(const std::vector<bead_ptr_t>& all,
               const std::vector<bead_ptr_t>& free,
               const std::vector<bead_group_ptr_t>& groups,
               const PairLists<Bead>& pairLists) = 0;
        
        
        /**
         * Returns the interaction energy of given bead with all other beads.
         * @param bead Bead
         * @param all All beads.
         * @param free Free beads.
         * @param groups bead groups.
         * @return Bonded and non-bonded interaction energy (potential energy).
         */
        virtual std::pair<energy_t, energy_t> 
        interact(const bead_ptr_t& bead,
                 const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups) = 0;
        
    };
}


#endif /* CG_FORCEFIELD_HPP */

