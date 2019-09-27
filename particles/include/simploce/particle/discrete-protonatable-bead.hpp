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
 * File:   discretely-protonatable-bead.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 5:33 PM
 */

#ifndef DISCRETE_PROTONATABLE_BEAD_HPP
#define DISCRETE_PROTONATABLE_BEAD_HPP

#include "bead.hpp"
#include "protonatable.hpp"
#include "ptypes.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Protonatable bead undergoing -discrete- or -discontinuous- changes in 
     * mass and charge values. Can bind more that one proton.
     */
    class DiscreteProtonatableBead: public Bead, Protonatable {
    public:
        
        virtual ~DiscreteProtonatableBead();
        
        charge_t charge() const override;
        
        mass_t mass() const override;
        
        void protonate() override;
        
        void deprotonate() override;
        
        bool isProtonated() const override;
                
        std::size_t protonationState() const override;
        
        void write(std::ostream& stream) const override;
        
        void writeState(std::ostream& stream) const override;
        
        void readState(std::istream& stream) override;        
        
    private:
        
        friend class CoarseGrained;
        
        static dprot_bead_ptr_t create(std::size_t id,
                                       std::size_t index, 
                                       const std::string& name,
                                       std::size_t numberOfBoundProtons,
                                       const spec_ptr_t& spec);
                
        DiscreteProtonatableBead(std::size_t id,
                                   std::size_t index, 
                                   const std::string &name,
                                   std::size_t numberOfBoundProtons,
                                   const spec_ptr_t &spec);
               
        std::size_t numberOfBoundProtons_;
    };
    
    /**
     * Writes a protonatable to an output stream.
     * @param stream Output stream.
     * @param pbead Protonatable.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, 
                               const DiscreteProtonatableBead& bead);
}

#endif /* PROTONATABLE_BEAD_HPP */

