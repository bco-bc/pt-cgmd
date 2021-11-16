/*
 * File:   pol-water-force-field.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 30, 2019, 2:02 PM
 */

#ifndef CG_POL_WATER_HPP
#define CG_POL_WATER_HPP

#include "cg-forcefield.hpp"
#include "sim-data.hpp"
#include "s-types.hpp"
#include <vector>

namespace simploce {
    
  /**
   * Force field for the polarizable coarse grained model according to 
   * Riniker and van Gunsteren, J. Chem. Phys. 134, 084119, 2011.
   */
    class CoarseGrainedPolarizableWater : public CoarseGrainedForceField {
    public:
        
        CoarseGrainedPolarizableWater(const spec_catalog_ptr_t& catalog,
                                      const bc_ptr_t& bc,
                                      const box_ptr_t& box,
                                      bool protonatable);
        
        std::pair<energy_t, energy_t>
        interact(const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups,
                 const PairLists<Bead>& pairLists) override;
        
        energy_t bonded(const std::vector<bead_ptr_t>& all,
                        const std::vector<bead_ptr_t>& free,
                        const std::vector<bead_group_ptr_t>& groups,
                        const PairLists<Bead>& pairLists) override; 
        
        std::pair<energy_t, energy_t>
        interact(const bead_ptr_t& bead,
                 const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups) override;
        
        std::string id() const override;
        
        /**
         * Includes C12, and C6 LJ parameters for the particle specifications 
         * pair (CW, CW).
         * @return Parameters.
         */
        std::pair<lj_params_t, el_params_t> parameters() const override;
        
        /**
         * Ideal distance between CW and DP.
         * @return Distance, in nm.
         */
        static length_t idealDistanceCWDP();
        
    private:
        
        spec_catalog_ptr_t catalog_;
        bc_ptr_t bc_;
        box_ptr_t box_;
    };
}

#endif /* POL_WATER_FORCE_FIELD_HPP */

