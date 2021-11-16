/*
 * File:   distance-lists.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 12:13 PM
 */

#ifndef DISTANCE_LISTS_HPP
#define DISTANCE_LISTS_HPP

#include "pair-list-generator.hpp"
#include "simploce/particle/bead.hpp"
#include "s-types.hpp"

namespace simploce {
    
    /**
     * Creates particle pair lists based on distances between particles.
     * @param P Particle type.
     */
    template <typename P>
    class DistanceLists;
    
    /**
     * Specialization for atoms
     */
    template <>
    class DistanceLists<Atom> : public pair_lists_generator<Atom> {
    public:
        
        DistanceLists(const box_ptr_t& box,
                      const bc_ptr_t& bc);
        
        PairLists<Atom>
        generate(const std::vector<atom_ptr_t>& all,
                 const std::vector<atom_ptr_t>& free,
                 const std::vector<atom_group_ptr_t>& groups) const override;
        
    private:
        
        box_ptr_t box_;
        bc_ptr_t bc_;
                
    };
    
    
    /**
     * Specialization for beads.
     */
    template <>
    class DistanceLists<Bead> : public pair_lists_generator<Bead> {
    public:
                        
        DistanceLists(const box_ptr_t& box,
                      const bc_ptr_t& bc);
        
        PairLists<Bead> 
        generate(const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups) const override;
        
    private:
        
        box_ptr_t box_;
        bc_ptr_t bc_;
                
    };
}


#endif /* DISTANCE_LISTS_HPP */

