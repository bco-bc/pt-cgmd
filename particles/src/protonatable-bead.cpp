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
 * File:   protonatable-bead.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 3:34 PM
 */

#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/particle/ptypes.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/pconf.hpp"
#include <stdexcept>

namespace simploce {
    
    void ProtonatableBead::protonate()
    {
        numberOfBoundProtons_ += 1;
    }
    
    void ProtonatableBead::deprotonate()
    {
        if ( !this->isProtonated() ) {
            throw std::domain_error(
                "Illegal attempt to deprotonate an already deprotonated bead."
            );
        }
        numberOfBoundProtons_ -= 1;              
    }
    
    bool ProtonatableBead::isProtonated() const
    {
        return numberOfBoundProtons_ > 0;
    }
    
    std::size_t ProtonatableBead::protonationState() const
    {
        return numberOfBoundProtons_;
    }
    
    ProtonatableBead::~ProtonatableBead()
    {        
    }
    
    void ProtonatableBead::write(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        Bead::write(stream);
        stream << space << this->protonationState();
    }
    
    void ProtonatableBead::writeState(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        Bead::writeState(stream);
        stream << space << this->protonationState();
    }
    
    ProtonatableBead::ProtonatableBead(std::size_t index, 
                                       const std::string& name,
                                       std::size_t numberOfBoundProtons,
                                       const bead_spec_ptr_t& spec) :
        Bead(index, name, spec), numberOfBoundProtons_(numberOfBoundProtons)
    {        
    }
    
    prot_bead_ptr_t ProtonatableBead::create(std::size_t index, 
                                             const std::string& name,
                                             std::size_t numberOfBoundProtons,
                                             const bead_spec_ptr_t& spec)
    {
        if ( !spec->isProtonatable() ) {
            throw std::domain_error(
                spec->name() + 
                ": this specification does not allow for (de)protonation."
            );
        }
        return prot_bead_ptr_t(new ProtonatableBead(index, 
                                                    name, 
                                                    numberOfBoundProtons, 
                                                    spec));
    }
    
    std::ostream& operator << (std::ostream& stream, const ProtonatableBead& pbead)
    {
        pbead.write(stream);
        return stream;
    }
    
}

