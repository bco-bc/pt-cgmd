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
         * Returns mass. If this specification is for a protonatable, then 
         * the value of the fully deprotonated state is returned.
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
         * This specification is for a protonatable particle, either with 
         * continuously or discretely varying protonation state) ?
         * @return Result.
         */
        bool isProtonatable() const;
        
        /**
         * This specification is for a protonatable particle with continuously varying 
         * protonation state.
         * @return Result.
         */
        bool isContinuousProtonatable() const;
        
        /**
         * Creates a new specification from a given specification, but with 
         * different name and charge.
         * @param spec Specification.
         * @return New specification.
         */
        static spec_ptr_t createFrom(const spec_ptr_t& spec, 
                                     const std::string& name,
                                     charge_t charge);
        
        /**
         * Creates a non-protonatable particle specification.
         * @param name Unique specification name.
         * @param charge Charge value.
         * @param mass Mass value.
         * @param radius Radius.
         * @return Specification.
         */
        static spec_ptr_t create(const std::string& name,
                                 charge_t charge,
                                 mass_t mass,
                                 radius_t radius);
        
        /**
         * Creates a protonatable particle specification.
         * @param name Unique specification name.
         * @param charge Charge value. This value must be that of the fully 
         * -deprotonated- state.
         * @param mass Mass value. This value must be that of the fully 
         * -deprotonated- state.
         * @param radius Radius.
         * @param pKa pKa value. 
         * @param continuous If true, then this specification is meant for
         * a protonatable with a continuously varying protonation state.
         * @return Specification.
         */
        static spec_ptr_t create(const std::string& name,
                                 charge_t charge,
                                 mass_t mass,
                                 radius_t radius,
                                 pKa_t pKa,
                                 bool continuous);
        
        /**
         * Writes this specification to an output stream.
         * @param stream Output steam.
         */
        void write(std::ostream& stream) const;
        
    private:
        
        ParticleSpec(const std::string& name,
                     charge_t charge,
                     mass_t mass,
                     radius_t radius,
                     pKa_t pKa,
                     bool protonatable,
                     bool continuous);
                
        std::string name_;
        charge_t charge_;
        mass_t mass_;
        radius_t radius_;
        pKa_t pKa_;        
        bool protonatable_;
        bool continuous_;        
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

