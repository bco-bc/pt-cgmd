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
 * File:   particle.hpp
 * Author: André H. Juffer, Biocenter Oulu
 *
 * Created on July 30, 2019, 4:05 PM
 */

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "ptypes.hpp"
#include <memory>
#include <string>

namespace simploce {
    
    /**
     * Recognizable unit composing the physical system at any time. Has location,
     * feels forces, and may be moving. Also, has charge and mass.
     */
    class Particle {
    public:
        
        // Noncopyable.
        Particle(const Particle&) = delete;
        Particle& operator = (const Particle&) = delete;
        
        virtual ~Particle();
        
        /**
         * Returns unique particle identifier.
         * @return Index, always >= 1.
         */
        std::size_t id() const;
        
        /**
         * Returns unique particle index.
         * @return Index, always >= 1.
         */
        std::size_t index() const;
        
        /**
         * Returns name. May not be unique.
         * @return Name.
         */
        std::string name() const;
        
        /**
         * Returns specification.
         * @return 
         */
        particle_spec_ptr_t spec() const;
        
        /**
         * Returns total charge value.
         * @return Value.
         */
        virtual charge_t charge() const;
        
        /**
         * Returns total mass value.
         * @return Value, always > 0.
         */
        virtual mass_t mass() const;
        
        /**
         * Returns position.
         * @return Position.
         */
        const position_t position() const;
        
        /**
         * Sets position.
         * @param r New position.
         */
        void position(const position_t& r);
        
        /**
         * Returns linear momentum.
         * @return Momentum.
         */
        const momentum_t momentum() const;
        
        /**
         * Set momentum.
         * @param p Momentum.
         */
        void momentum(const momentum_t& p);
        
        /**
         * Returns force acting on this particle.
         * @return Force.
         */
        const force_t force() const;
        
        /**
         * Sets force acting on this particle.
         * @param f Force.
         */
        void force(const force_t& f);
        
        /**
         * Resets force acting on this particle to zero.
         */
        void resetForce();
        
        /**
         * Writes this particle to an output stream.
         * @param stream Output stream.
         */
        virtual void write(std::ostream& stream) const;
        
        /**
         * Writes state to an output stream.
         * @param stream Output stream.
         */
        virtual void writeState(std::ostream& stream) const;
        
    protected:
        
        /**
         * Constructor
         * @param id Unique particle index.
         * @param name Particle name. Does not need to be unique.
         * @param spec Particle specification.
         */
        Particle(std::size_t index, 
                 const std::string& name, 
                 const particle_spec_ptr_t& spec);        
                        
    private:
        
        template <typename P, typename PG>
        friend class ParticleModel;

        template <typename T, typename S>
        friend class ProtonationSite;
        
        void reset_(const particle_spec_ptr_t& spec);
        
        //virtual bool isProtonatable_() const;
        
        //virtual std::size_t protonationState_() const;
        
        std::size_t index_;
        std::string name_;
        particle_spec_ptr_t spec_;
        position_t r_;
        momentum_t p_;
        force_t f_;
        
    };
    
    /**
     * Writes particle to output stream.
     * @param stream Output stream.
     * @param particle Particle.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const Particle& particle);
}



#endif /* PARTICLE_HPP */

