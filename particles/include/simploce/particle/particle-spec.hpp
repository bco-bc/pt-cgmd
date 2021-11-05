/*
 * File:   particle-spec.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 4:24 PM
 */

#ifndef PARTICLE_SPEC_HPP
#define PARTICLE_SPEC_HPP

#include "p-types.hpp"
#include <string>
#include <iostream>

namespace simploce {
    
    /**
     * Defines certain static properties of particles. Multiple particles may 
     * hold or refer to the same specification.
     */
    class ParticleSpec {
    public:

        /**
         * Creates a non-protonatable particle specification. All arguments are required.
         * @param name Unique specification name.
         * @param charge Charge value.
         * @param mass Mass value. Must be a non-negative number.
         * @param radius Radius. Must be a non-negative number.
         * @param description Short description.
         * @return Specification.
         */
        static spec_ptr_t create(const std::string& name,
                                 const charge_t& charge,
                                 const mass_t& mass,
                                 const radius_t& radius,
                                 const std::string& description);

        /**
         * Creates a protonatable particle specification.
         * @param name Unique specification name.
         * @param charge Charge value. Must represent the -deprotonated- state.
         * @param mass Mass value. Must be a non-negative number and must represent the -deprotonated- state.
         * @param radius Radius.  Must be a non-negative number.
         * @param pKa pKa value. Typically in the range [0.0, 14.0].
         * @param description Short description. Must be provided.
         * @return Specification.
         */
        static spec_ptr_t create(const std::string& name,
                                 const charge_t& charge,
                                 const mass_t& mass,
                                 const radius_t& radius,
                                 const pKa_t& pKa,
                                 const std::string& description);

        /**
         * Creates a new specification from a given specification, but with
         * different unique name and charge.
         * @param spec Specification.
         * @param name New name.
         * @param charge New Charge value.
         * @param description Short description.
         * @return New specification.
         */
        static spec_ptr_t createFrom(const spec_ptr_t& spec,
                                     const std::string& name,
                                     charge_t charge,
                                     const std::string& description);

        // Noncopyable.
        ParticleSpec(const ParticleSpec&) = delete;
        ParticleSpec& operator = (const ParticleSpec&) = delete;

        // Compares unique names.
        bool operator == (const ParticleSpec& particleSpec);

        /**
         * Returns unique identifying name.
         * @return Name.
         */
        std::string name() const;
        
        /**
         * Returns charge. If representing a protonatable specification, then the value
         * of the deprotonated state is returned.
         * @return Value.
         * @see #isProtonatable()
         */
        charge_t charge() const;
        
        /**
         * Returns mass. If representing a protonatable specification, then the value
         * of the deprotonated state is returned.
         * @return Value.
         * @see #isProtonatable()
         */
        mass_t mass() const;
        
        /**
         * Returns pKa value.
         * @return pKa.
         */
        pKa_t pKa() const;
        
        /**
         * Returns radius.
         * @return Radius.
         */
        radius_t radius() const;
        
        /**
         * Is this specification for a protonatable particle.
         * @return Result.
         */
        bool isProtonatable() const;
        
        /**
         * Is this a specification for an ion.
         * @return Result.
         */
        bool isIon();

        /**
         * Returns description.
         * @return Description.
         */
        std::string description() const;
        
        /**
         * Writes this specification to an output stream.
         * @param stream Output steam.
         */
        void write(std::ostream& stream) const;
        
    private:

        static spec_ptr_t create_(const std::string& name,
                                  const charge_t& charge,
                                  const mass_t& mass,
                                  const radius_t& radius,
                                  const pKa_t& pKa,
                                  bool protonatable,
                                  const std::string& description);

        ParticleSpec(std::string name,
                     charge_t charge,
                     mass_t mass,
                     radius_t radius,
                     pKa_t pKa,
                     bool protonatable,
                     std::string description);

        std::string name_;
        charge_t charge_;
        mass_t mass_;
        radius_t radius_;
        pKa_t pKa_;        
        bool protonatable_;
        std::string description_;
    };

    /**
     * Write particle specification to output stream.
     * @param stream Output stream.
     * @param spec Particle specification.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const ParticleSpec& spec);
}


#endif /* PARTICLE_SPEC_HPP */

