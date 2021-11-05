/*
 * File:   atom.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 4:46 PM
 */

#ifndef ATOM_HPP
#define ATOM_HPP

#include "particle.hpp"
#include "p-types.hpp"
#include <string>

namespace simploce {
    
    /**
     * The smallest constituent unit of ordinary matter that has the 
     * properties of a chemical element. An atom cannot freely be created.
     * Each atom belongs to an Atomistic
     * @see Atomistic
     */
    class Atom: public Particle {
    public:

        /**
         * Creates a new atom. All arguments are required.
         * @param id Unique atom identifier.
         * @param index Sequential index.
         * @param name Atom name.
         * @param spec Atom specification.
         * @return
         */
        static atom_ptr_t create(const id_t& id,
                                 std::size_t index,
                                 const std::string &name,
                                 const spec_ptr_t &spec);

        /**
         * Constructor. All arguments are required.
         * @param id Unique atom identifier.
         * @param index Sequential index.
         * @param name Atom name.
         * @param spec Atom specification.
         */
        Atom(const id_t& id,
             std::size_t index,
             const std::string &name,
             const spec_ptr_t &spec);

        // Noncopyable.
        Atom(const Atom &) = delete;
        Atom &operator=(const Atom &) = delete;

        virtual ~Atom();

    };

}

#endif /* ATOM_HPP */

