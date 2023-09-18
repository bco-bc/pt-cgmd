/*
 * File:   bead.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 4:44 PM
 */

#ifndef BEAD_HPP
#define BEAD_HPP

#include "particle.hpp"
#include "p-types.hpp"

namespace simploce {
    
    /**
     * Coarse grained particle. In proteins, it represents typically 4 to 6 atoms,
     * but for water it may represent several water molecules.
     * @see CoarseGrained
     */
    class Bead : public Particle {
    public:

        /**
         * Creates a new bead. All arguments are required.
         * @param id Unique bead identifier.
         * @param index Sequential index.
         * @param name Bead name.
         * @param spec Bead specification. Must -not- refer to a protonatable.
         * @return Bead.
         */
        static bead_ptr_t create(const id_t& id,
                                 std::size_t index,
                                 const std::string& name,
                                 const spec_ptr_t& spec);

        // Noncopyable.
        Bead(const Bead&) = delete;
        Bead& operator = (const Bead&) = delete;
        
        ~Bead() override;

    protected:

        /**
         * Constructor. All arguments are required.
         * @param id Unique bead identifier.
         * @param index Sequential index.
         * @param name Bead name.
         * @param spec Bead specification. Must -not- refer to a protonatable.
         */
        Bead(const id_t& id,
             std::size_t index,
             const std::string& name,
             const spec_ptr_t& spec);
    };
}


#endif /* BEAD_HPP */

