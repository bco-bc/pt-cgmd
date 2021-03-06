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
    
    Particle::Particle(std::size_t id,
                       std::size_t index, 
                       const std::string& name, 
                       const spec_ptr_t& spec) :
        id_{id}, index_{index}, name_{name}, spec_{spec}, r_{}, p_{}, v_{}, 
        f_{}, pf_{}
    {
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
    
    std::size_t 
    Particle::id() const
    {
        return id_;
    }
    
    std::size_t 
    Particle::index() const
    {
        return index_;
    }
    
    std::string 
    Particle::name() const
    {
        return name_;
    }
    
    spec_ptr_t 
    Particle::spec() const
    {
        return spec_;
    }
    
    charge_t 
    Particle::charge() const
    {
        return spec_->charge();
    }
    
    mass_t 
    Particle::mass() const
    {
        return spec_->mass();
    }
    
    const 
    position_t Particle::position() const
    {
        return r_;
    }
    
    void 
    Particle::position(const position_t& r) 
    {
        r_ = r; 
    }
    
    const 
    momentum_t Particle::momentum() const
    { 
        real_t ma = this->mass()();
        velocity_t vv = ma * v_;
        return momentum_t{vv.toArray()};
    }
    
    void 
    Particle::momentum(const momentum_t& p) 
    { 
        real_t ma = this->mass()();
        for (std::size_t k = 0; k != 3; ++k) {
            v_[k] = p[k] / ma;
        }
    }
    
    velocity_t 
    Particle::velocity() const
    {
        return v_;
    }
    
    void 
    Particle::velocity(const velocity_t& v)
    {
        v_ = v;
    }
    
    const 
    force_t Particle::force() const 
    {
        return f_; 
    }
    
    void 
    Particle::force(const force_t& f) 
    { 
        f_ = f;
    }
    
    void 
    Particle::resetForce() 
    { 
        f_ = force_t{}; 
    }
    
    void 
    Particle::write(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        stream << std::setw(10) << this->name();
        stream << std::setw(10) << this->spec()->name();
        stream << space << std::setw(10) << this->id();
        stream << space << this->position();
        stream << space << this->velocity();
    }
    
    void 
    Particle::writeState(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        stream << space << this->position();
        stream << space << this->velocity();
    }
    
    void 
    Particle::readState(std::istream& stream)
    {
        real_t x, y, z, vx, vy, vz;
        stream >> x >> y >> z >> vx >> vy >> vz;
        position_t r{x, y, z};
        velocity_t v{vx, vy, vz};
        this->position(r);
        this->velocity(v);
    }
    
    void 
    Particle::reset_(const spec_ptr_t &spec)
    {
        if ( !spec ) {
            throw std::domain_error("A particle specification must be provided.");
        }
        spec_ = spec;
    }
    
    std::ostream& 
    operator << (std::ostream& stream, 
                 const Particle& particle)
    {
        particle.write(stream);
        return stream;
    }
    
}
