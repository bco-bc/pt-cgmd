/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* 
 * File:   particle-spec.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 4:24 PM
 */

#ifndef PARTICLE_SPEC_HPP
#define PARTICLE_SPEC_HPP

#include "ptypes.hpp"
#include <string>
#include <iostream>

namespace simploce {
    
    /**
     * Defines certain static properties of particles. Multiple particles may 
     * hold the same specification.
     */
    class ParticleSpec {
    public:
        
        /**
         * Returns unique identifying name.
         * @return Name.
         */
        std::string name() const;
        
        /**
         * Returns charge. If protonatable, then the value of the 
         * deprotonated state is returned.
         * @return Value.
         * @see #isProtonatable()
         */
        charge_t charge() const;
        
        /**
         * Returns mass. If protonatable, then the value of the 
         * deprotonated state is returned.
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
         * Particle size.
         * @return Radius.
         */
        radius_t radius() const;
        
        /**
         * This specification allows for binding and releasing protons?
         * @return Result.
         */
        bool isProtonatable() const;
        
        /**
         * Creates a new specification from a given specification, but with different 
         * name and charge.
         * @param spec Specification.
         * @return New specification.
         */
        static particle_spec_ptr_t createFrom(const particle_spec_ptr_t& spec, 
                                              const std::string& name,
                                              charge_t charge);
        
        /**
         * Creates a non-protonatable bead specification.
         * @param name Unique specification name.
         * @param charge Charge value.
         * @param mass Mass value.
         * @param radius Radius.
         * @return Bead Specification.
         */
        static bead_spec_ptr_t createForBead(const std::string& name,
                                             charge_t charge,
                                             mass_t mass,
                                             radius_t radius);
        
        /**
         * Creates a protonatable bead specification.
         * @param name Unique specification name.
         * @param charge Charge value. This value must be that of the fully 
         * -deprotonated- state.
         * @param mass Mass value. This value must be that of the fully 
         * -deprotonated- state.
         * @param radius Radius.
         * @param pKa pKa value. 
         * @return Protonatable bead Specification.
         */
        static prot_bead_spec_ptr_t createForProtonatableBead(const std::string& name,
                                                              charge_t charge,
                                                              mass_t mass,
                                                              radius_t radius,
                                                              pKa_t pKa);
        
        /**
         * Creates an atom specification. Always not protonatable.
         * @param name Unique specification name.
         * @param charge Charge value.
         * @param mass Mass value.
         * @param radius Radius.
         * @return Atom specification.
         */
        static atom_spec_ptr_t createForAtom(const std::string& name,
                                             charge_t charge,
                                             mass_t mass,
                                             radius_t radius);
        
        /**
         * Writes this specification to an output stream.
         * @param stream Output steam.
         */
        void writeTo(std::ostream& stream) const;
        
    private:
        
        ParticleSpec(const std::string& name,
                     charge_t charge,
                     mass_t mass,
                     radius_t radius,
                     pKa_t pKa,
                     bool protonatable,
                     const std::string& type);
                
        std::string name_;
        charge_t charge_;
        mass_t mass_;
        radius_t radius_;
        pKa_t pKa_;        
        bool protonatable_;
        
        // One of "atom", "bead" or "prot_bead".
        std::string type_;
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

