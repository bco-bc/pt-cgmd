/*
 * The MIT License
 *
 * Copyright 2019 juffer.
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
 * File:   protonating-bead.hpp
 * Author: juffer
 *
 * Created on 25 September 2019, 20:25
 */

#ifndef PROTONATING_BEAD_HPP
#define PROTONATING_BEAD_HPP

#include "bead.hpp"
#include "protonatable-bead.hpp"
#include "ptypes.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Protonatable undergoing in - continuous- changes in mass and charge values. 
     * This is due to involvement in proton transfer processes, the state of 
     * which is measured from a current or flow of mass and charge to and from 
     * the bead. If the state is denoted as x, then the associated 
     * current I = dx/dt, where t is time. Here, x is in [0,1], where a value 
     * of 1(0) represents a fully (de)protonated state.
     */
    class ProtonatingBead : public Bead, Protonatable  {
    public:
        
        virtual ~ProtonatingBead();
        
        /**
         * Returns state.
         * @return State.
         */
        real_t state() const;
        
        /**
         * Sets state.
         * @param x State. Must be in [0,1].
         */
        void state(real_t x);
        
        /**
         * Returns current.
         * @return Current.
         */
        real_t current() const;
        
        /**
         * Sets current.
         * @param I Current.
         */
        void current(real_t I);
        
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
        
        static cprot_bead_ptr_t create(std::size_t id,
                                       std::size_t index, 
                                       const std::string& name,
                                       std::size_t protonationState,
                                       const bead_spec_ptr_t& spec);
        
        ProtonatingBead(std::size_t id,
                        std::size_t index, 
                        const std::string &name,
                        std::size_t protonationState,
                        const bead_spec_ptr_t &spec);
        
        // Represents state of transfer, relative to a fully deprotonated state.
        real_t x_;
        real_t I_;  // Current, I = dx/dt.
        std::size_t protonationState_;
        
    };
    
    std::ostream& operator << (std::ostream& stream, const ProtonatingBead& bead);
}

#endif /* PROTONATING_BEAD_HPP */

