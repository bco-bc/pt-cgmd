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
 * File:   discretely-protonatable-bead.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 3:34 PM
 */

#include "simploce/particle/discrete-protonatable-bead.hpp"
#include "simploce/particle/ptypes.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/pconf.hpp"
#include <stdexcept>

namespace simploce {
    
    DiscreteProtonatableBead::~DiscreteProtonatableBead()
    {        
    }
    
    charge_t DiscreteProtonatableBead::charge() const 
    {
        // Fully deprotonated state.
        charge_t charge = Particle::charge();
        
        // Correct for bound proton.
        charge += this->protonationState() * conf::CHARGE_PROTON;
        return charge;
    }
    
    mass_t DiscreteProtonatableBead::mass() const
    {
        // Fully deprotonated state.
        mass_t mass = Particle::mass();
        
        // Correct for bound proton.
        mass += this->protonationState() * conf::MASS_PROTON;
        
        return mass;
    }
    
    void DiscreteProtonatableBead::protonate()
    {
        numberOfBoundProtons_ += 1;
    }
    
    void DiscreteProtonatableBead::deprotonate()
    {
        if ( !this->isProtonated() ) {
            throw std::domain_error(
                "Illegal attempt to deprotonate an already fully deprotonated bead."
            );
        }
        numberOfBoundProtons_ -= 1;              
    }
    
    bool DiscreteProtonatableBead::isProtonated() const
    {
        return numberOfBoundProtons_ > 0;
    }
    
    std::size_t DiscreteProtonatableBead::protonationState() const
    {
        return numberOfBoundProtons_;
    }
    
    void DiscreteProtonatableBead::write(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        Particle::write(stream);
        stream << space << conf::DISCRETELY_PROTONATABLE << space << numberOfBoundProtons_;
    }
    
    void DiscreteProtonatableBead::writeState(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        Particle::writeState(stream);
        stream << space << this->protonationState();
    }
    
    void DiscreteProtonatableBead::readState(std::istream& stream)
    {
        Particle::readState(stream);
        stream >> numberOfBoundProtons_;        
    }
    
    dprot_bead_ptr_t DiscreteProtonatableBead::create(std::size_t id,
                                              std::size_t index, 
                                              const std::string& name,
                                              std::size_t numberOfBoundProtons,
                                              const spec_ptr_t& spec)
    {
        if ( !spec->isProtonatable() ) {
            throw std::domain_error(
                spec->name() + 
                ": this specification does not allow for (de)protonation."
            );
        }
        return dprot_bead_ptr_t(new DiscreteProtonatableBead(id,
                                                     index, 
                                                     name, 
                                                     numberOfBoundProtons, 
                                                     spec));
    }
    
    DiscreteProtonatableBead::DiscreteProtonatableBead(std::size_t id,
                                                       std::size_t index, 
                                                       const std::string& name,
                                                       std::size_t numberOfBoundProtons,
                                                       const spec_ptr_t& spec) :
        Bead(id, index, name, spec), numberOfBoundProtons_(numberOfBoundProtons)
    {        
    }
    
    std::ostream& operator << (std::ostream& stream, 
                               const DiscreteProtonatableBead& bead)
    {
        bead.write(stream);
        return stream;
    }
    
}

