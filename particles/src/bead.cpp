/*
 * File:   bead.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 2:00 PM
 */

#include "simploce/particle/bead.hpp"

namespace simploce {
    
    Bead::Bead(const id_t& id,
               std::size_t index, 
               const std::string& name, 
               const spec_ptr_t& spec) :
        Particle(id, index, name, spec) {
    }
        
    Bead::~Bead() = default;
    
    bead_ptr_t
    Bead::create(const id_t& id,
                 std::size_t index,
                 const std::string& name,
                 const spec_ptr_t& spec) {
        return std::make_shared<Bead>(id, index, name, spec);
    }
        
}

