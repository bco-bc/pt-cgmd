/*
 * File:   coarse-grained.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:34 PM
 */

#ifndef COARSE_GRAINED_HPP
#define COARSE_GRAINED_HPP

#include "particle-system.hpp"
#include "p-types.hpp"

namespace simploce {
    
    /**
     * A particle system composed of regular non-protonatable beads.
     */
    class CoarseGrained : public ParticleSystem<Bead, bead_group_t> {
    public:

        /**
         * Reads a coarse grained particle system from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specifications catalog.
         * @return Coarse grained particle system.
         */
        static cg_mod_ptr_t obtainFrom(std::istream& stream,
                                       const spec_catalog_ptr_t& catalog);


        /**
         * Constructor. Creates empty coarse grained particle system.
         */
        CoarseGrained();
        
        /**
         * Adds a new bead to this physical system. All arguments are required.
         * @param name Name. Not required to be unique.
         * @param spec Bead specification.
         * @return Newly created bead.
         */
        bead_ptr_t addBead(const std::string& name,
                           const spec_ptr_t& spec);

        /**
         * Adds a bead group with bonds to this physical system. All arguments are required.
         * @param beads Beads forming a bead group. Each of these beads must already
         * be present in this particle model.
         * @param bonds Bonds between beads, given as pairs of bead identifiers.
         * @return Bead group added.
         */
        bead_group_ptr_t addBeadGroup(const std::vector<bead_ptr_t>& beads, 
                                      const std::vector<id_pair_t>& bonds);
        
        /**
         * Returns number of beads.
         * @return Number.
         */
        std::size_t numberOfBeads() const;

    private:

        bead_ptr_t createParticle_(const id_t& id,
                                   int index,
                                   const std::string& name,
                                   const spec_ptr_t& spec) override;

        friend class ParticleModelFactory;

    };

    /**
     * Writes coarse grained particle system to output stream.
     * @param stream Output stream.
     * @param cg Coarse grained particle system.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const CoarseGrained& cg);
    
}

#endif /* COARSE_GRAINED_HPP */

