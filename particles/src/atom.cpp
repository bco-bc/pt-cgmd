/*
 * File:   atom.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 2:58 PM
 */

#include "simploce/particle/atom.hpp"

namespace simploce {
    
    Atom::Atom(const id_t& id,
               std::size_t index, 
               const std::string& name, 
               const spec_ptr_t& spec) :
        Particle(id, index, name, spec) {
    }
    
    Atom::~Atom() = default;
    
    atom_ptr_t 
    Atom::create(const id_t& id,
                 std::size_t index, 
                 const std::string &name, 
                 const spec_ptr_t &spec) {
        return atom_ptr_t(new Atom(id, index, name, spec));
    }
}