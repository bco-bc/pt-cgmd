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
 * File:   particle-spec.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 2:01 PM
 */

#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/pconf.hpp"
#include <memory>
#include <stdexcept>

namespace simploce {
    
    std::string 
    ParticleSpec::name() const
    {
        return name_;
    }
    
    charge_t 
    ParticleSpec::charge() const
    {
        return charge_;
    }
    
    mass_t 
    ParticleSpec::mass() const
    {
        return mass_;
    }
    
    radius_t 
    ParticleSpec::radius() const
    {
        return radius_;
    }
    
    pKa_t 
    ParticleSpec::pKa() const
    {
        return pKa_;
    }
    
    bool 
    ParticleSpec::isProtonatable() const
    {
        return protonatable_;
    }
    
    bool 
    ParticleSpec::isContinuousProtonatable() const
    {
        return continuous_;
    }
    
    spec_ptr_t 
    ParticleSpec::createFrom(const spec_ptr_t& spec,
                             const std::string& name,
                             charge_t charge)
    {
        return spec_ptr_t(new ParticleSpec(
            name, charge, spec->mass(), spec->radius(), spec->pKa(), 
            spec->isProtonatable(), spec->isContinuousProtonatable()
        ));
    }
    
    spec_ptr_t 
    ParticleSpec::create(const std::string& name, 
                         charge_t charge, 
                         mass_t mass,
                         radius_t radius)
    {
        return spec_ptr_t(
            new ParticleSpec(name, charge, mass, radius, 0.0, false, false)
        );
    }
    
    spec_ptr_t 
    ParticleSpec::create(const std::string& name,
                         charge_t charge,
                         mass_t mass,
                         radius_t radius,
                         pKa_t pKa,
                         bool continuous)
    {
        return spec_ptr_t(
            new ParticleSpec(name, charge, mass, radius, pKa, true, continuous)
        );
    }
        
    void 
    ParticleSpec::write(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        const auto width = conf::NAME_WIDTH;
        
        stream << std::setw(width) << protonatable_;
        stream << std::setw(width) << name_;
        stream << space << mass_;
        stream << space << charge_;
        stream << space << radius_;
        if ( this->isProtonatable() ) {
            stream << space << pKa_ << space << continuous_;
        }        
    }
    
    ParticleSpec::ParticleSpec(const std::string& name,
                               charge_t charge,
                               mass_t mass,
                               radius_t radius,
                               pKa_t pKa, 
                               bool protonatable,
                               bool continuous) :
        name_{name}, charge_{charge}, mass_{mass}, radius_{radius}, pKa_{pKa}, 
        protonatable_{protonatable}, continuous_{continuous}
    {        
        if ( name_.empty() ) {
            throw std::domain_error("A particle specification name must be provided.");
        }
        if ( mass_ <= 0.0) {
            throw std::domain_error("A mass must be a positive number.");
        }        
    }
        
    std::ostream& 
    operator << (std::ostream& stream, 
                 const ParticleSpec& spec)
    {
        spec.write(stream);
        return stream;
    }    
}
