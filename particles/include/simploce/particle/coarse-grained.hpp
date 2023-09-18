/*
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
    class CoarseGrained : public ParticleSystem {
    public:

        /**
         * Reads a coarse grained particle system from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specifications catalog.
         * @return Coarse grained particle system.
         */
        static cg_sys_ptr_t parseIt(std::istream& stream,
                                    const spec_catalog_ptr_t& catalog);


        /**
         * Constructor. Creates empty coarse grained particle system.
         */
        CoarseGrained();

        ~CoarseGrained() noexcept override;
        
        /**
         * Adds a new bead to this physical system. All arguments are required.
         * @param name Name. Not required to be unique.
         * @param spec Bead specification.
         * @return Newly created bead.
         */
        p_ptr_t addBead(const std::string& name,
                        const spec_ptr_t& spec);

        /**
         * Returns number of beads.
         * @return Number.
         */
        std::size_t numberOfBeads() const;

    private:

        p_ptr_t createParticle(const id_t& id,
                               int index,
                               const std::string& name,
                               const spec_ptr_t& spec) override;

        friend class ParticleSystemFactory;

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

