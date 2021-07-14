/*
 * File:   atom.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 4:46 PM
 */

#ifndef ATOM_HPP
#define ATOM_HPP

#include "particle.hpp"
#include "ptypes.hpp"
#include <string>

namespace simploce {
    
    /**
     * The smallest constituent unit of ordinary matter that has the 
     * properties of a chemical element.
     */
    class Atom: public Particle {
    public:
        
        // Noncopyable.
        Atom(const Atom&) = delete;
        Atom& operator = (const Atom&) = delete;
        
        ~Atom();
                
    private:
        
        friend class Atomistic;
        
        static atom_ptr_t create(std::size_t id,
                                 std::size_t index, 
                                 const std::string &name, 
                                 const spec_ptr_t &spec);
        
        Atom(std::size_t id,
             std::size_t index, 
             const std::string &name,
             const spec_ptr_t &spec);        
        
    };
}

#endif /* ATOM_HPP */

