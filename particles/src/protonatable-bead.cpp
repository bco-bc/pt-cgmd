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
#include <stdexcept>

namespace simploce {
    
    prot_bead_ptr_t ProtonatableBead::create(int id, 
                                             const std::string& name,
                                             int protonationState,
                                             const bead_spec_ptr_t& spec)
    {
        if ( !spec->isProtonatable() ) {
            throw std::domain_error(
                spec->name() + 
                ": this specification does not allow for (de)protonation."
            );
        }
        return prot_bead_ptr_t(new ProtonatableBead(id, name, protonationState, spec));
    }
    
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
    
    ProtonatableBead::ProtonatableBead(int id, 
                                       const std::string& name,
                                       int protonationState,
                                       const bead_spec_ptr_t& spec) :
        Bead(id, name, spec), numberOfBoundProtons_(protonationState)
    {        
    }
    
    ProtonatableBead::~ProtonatableBead()
    {        
    }
    
    bool ProtonatableBead::isProtonatable_() const
    {
        return true;
    }
    
    std::size_t ProtonatableBead::protonationState_() const
    {
        return this->protonationState();
    }
    
}

