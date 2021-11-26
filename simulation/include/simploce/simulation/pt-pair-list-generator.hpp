/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 25, 2019, 1:05 PM
 */

#ifndef PT_PAIR_LIST_GENERATOR_HPP
#define PT_PAIR_LIST_GENERATOR_HPP

#include "s-types.hpp"
#include <memory>
#include <utility>

namespace simploce {
    
    /**
     * Finds pairs of protonatable beads that may be involved in transferring
     * protons.
     */
    class ProtonTransferPairListGenerator {
    public:
                
        /**
         * Protonatable beads pair typeName.
         */
        using prot_pair_t = std::pair<prot_bead_ptr_t, prot_bead_ptr_t>;
        
        /**
         * Protonatable bead pair list typeName.
         */
        using prot_pair_list_t = std::vector<prot_pair_t>;
        
        /**
         * Constructor 
         * @param bc Boundary condition.
         */
        ProtonTransferPairListGenerator(const bc_ptr_t& bc);
        
        /**
         * Generates protonatable bead pair list.
         * @param cg Coarse grained particle model.
         * @return List of pairs of beads possibly involved in proton transfer.
         */
        prot_pair_list_t generate(const prot_cg_sys_ptr_t& cg) const;
        
    private:
        
        length_t rcutoffPT_;
        bc_ptr_t bc_;
        
    };
    
     
}

#endif /* PT_PAIR_LIST_GENERATOR_HPP */

