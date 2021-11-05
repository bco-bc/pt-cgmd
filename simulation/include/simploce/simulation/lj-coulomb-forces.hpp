/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 17, 2019, 3:06 PM
 */

#ifndef LJ_COULOMB_HPP
#define LJ_COULOMB_HPP

#include "cg-forcefield.hpp"
#include "s-types.hpp"
#include "simploce/util/map2.hpp"
#include <string>
#include <map>
#include <vector>

namespace simploce {
    
    /**
     * Calculates LJ and Coulomb interaction.
     * @param P Particle type.
     */
    template <typename P>
    class LJCoulombForces;
    
    /**
     * Specialization for beads.
     */
    template <>
    class LJCoulombForces<Bead> : public CoarseGrainedForceField {
    public:
        
        LJCoulombForces(const lj_params_t& ljParams, 
                        const el_params_t& elParams,
                        const bc_ptr_t& bc,
                        const box_ptr_t& box);
        
        std::pair<energy_t, energy_t> 
        interact(const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups,
                 const PairLists<Bead>& pairLists) override;
        
        std::pair<energy_t, energy_t>
        interact(const bead_ptr_t& bead,
                 const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups) override;
        
        /**
         * There are non bonded interactions.
         * @return 0.0.
         */
        energy_t bonded(const std::vector<bead_ptr_t>& all,
                        const std::vector<bead_ptr_t>& free,
                        const std::vector<bead_group_ptr_t>& groups,
                        const PairLists<Bead>& pairLists) override;
        
        std::string id() const override;
        
        std::pair<lj_params_t, el_params_t> parameters() const override;
        
    private:
        
        lj_params_t ljParams_;
        el_params_t elParams_;
        bc_ptr_t bc_;
        box_ptr_t box_;
    };
}

#endif /* LJ_COULOMB_HPP */

