/*
 * File:   atomistic.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:19 PM
 */

#ifndef ATOMISTIC_HPP
#define ATOMISTIC_HPP

#include "particle-system.hpp"
#include "p-types.hpp"

namespace simploce {

    /**
     * A physical system composed of atoms and atom-like particles, such as ions.
     */
    class Atomistic : public ParticleSystem<Atom, atom_group_t> {
    public:

        /**
         * Reads an atomistic particle system from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specifications catalog.
         * @return Atomistic particle model system.
         */
        static at_mod_ptr_t obtainFrom(std::istream& stream,
                                       const spec_catalog_ptr_t& catalog);

        /**
         * Constructor. Creates empty atomistic particle model.
         */
        Atomistic();
        
        /**
         * Adds new atom to this physical system. All arguments are required.
         * @param name Name. Not required to be unique.
         * @param spec Specification.
         * @return Newly created atom.
         */
        atom_ptr_t addAtom(const std::string& name,
                           const spec_ptr_t& spec);
        
        /**
         * Returns number of atoms.
         * @return Number.
         */
        std::size_t numberOfAtoms() const;
        
    private:

        atom_ptr_t createParticle_(const id_t& id,
                                   int index,
                                   const std::string& name,
                                   const spec_ptr_t& spec) override;

        friend class ParticleModelFactory;
        
    };

    /**
     * Writes atomistic particle system to an output stream.
     * @param stream Output stream.
     * @param atomistic Atomistic particle system.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const Atomistic& atomistic);
}

#endif /* ATOMISTIC_HPP */

