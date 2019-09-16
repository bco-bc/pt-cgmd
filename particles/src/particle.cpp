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

#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <iomanip>

namespace simploce {        
    
    Particle::Particle(int id, const std::string& name, const particle_spec_ptr_t& spec) :
        id_{id}, name_{name}, spec_{spec}, r_{}, p_{}, f_{}
    {
        if ( id_ < 0 ) {
            throw std::domain_error(
                boost::lexical_cast<std::string, int>(id_) + ": Illegal particle identifier."
            );
        }
        if ( name_.empty() ) {
            throw std::domain_error("A particle name must be provided.");
        }
        if ( !spec_ ) {
            throw std::domain_error("A particle specification must be provided.");
        }
    }
    
    Particle::~Particle()
    {        
    }
    
    int Particle::id() const
    {
        return id_;
    }
    
    std::string Particle::name() const
    {
        return name_;
    }
    
    particle_spec_ptr_t Particle::spec() const
    {
        return spec_;
    }
    
    charge_t Particle::charge() const
    {
        return spec_->charge();
    }
    
    mass_t Particle::mass() const
    {
        return spec_->mass();
    }
    
    const position_t Particle::position() const
    {
        return r_;
    }
    
    void Particle::position(const position_t& r) 
    {
        r_ = r; 
    }
    
    const momentum_t Particle::momentum() const
    { 
        return p_; 
    }
    
    void Particle::momentum(const momentum_t& p) 
    { 
        p_ = p; 
    }
    
    const force_t Particle::force() const 
    {
        return f_; 
    }
    
    void Particle::force(const force_t& f) 
    { 
        f_ = f;
    }
    
    void Particle::resetForce() 
    { 
        f_ = force_t{}; 
    }
    
    void Particle::reset_(const particle_spec_ptr_t &spec)
    {
        if ( !spec ) {
            throw std::domain_error("A particle specification must be provided.");
        }
        spec_ = spec;
    }
    
    bool Particle::isProtonatable_() const
    {
        return false;
    }
    
    std::size_t Particle::protonationState_() const
    {
        return 0;
    }
    
    std::ostream& operator << (std::ostream& stream, const Particle& particle)
    {
        const char space = ' ';
        stream << std::setw(10) << particle.id();
        stream << space << std::setw(10) << particle.name();
        stream << space << std::setw(10) << particle.spec()->name();
        stream << space << particle.position();
        return stream;
    }
    
}
