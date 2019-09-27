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

#include "simploce/particle/continuous-protonatable-bead.hpp"
#include "simploce/particle/pconf.hpp"
#include <stdexcept>
#include <iomanip>
#include <cassert>

namespace simploce {
    
    ContinuousProtonatableBead::~ContinuousProtonatableBead()
    {        
    }
    
    real_t ContinuousProtonatableBead::state() const
    {
        return x_;
    }
    
    void ContinuousProtonatableBead::state(real_t x)
    {
        assert(x >= 0.0 && x <= 1.0);
        x_ = x;
    }
    
    real_t ContinuousProtonatableBead::current() const
    {
        return I_;
    }
    
    void ContinuousProtonatableBead::current(real_t I)
    {
        I_ = I;
    }
    
    charge_t ContinuousProtonatableBead::charge() const
    {
        // Fully deprotonated state.
        charge_t charge = Particle::charge();
        
        // From state of transfer.
        charge += x_ * conf::CHARGE_PROTON;
        
        return charge;
    }
    
    mass_t ContinuousProtonatableBead::mass() const
    {
        // Fully deprotonated state.
        mass_t mass = Particle::mass();
        
        // From state of proton transfer process.
        mass += x_ * conf::MASS_PROTON;
        
        return mass;
    }
    
    void ContinuousProtonatableBead::protonate()
    {
        protonationState_ = 1;
    }
    
    void ContinuousProtonatableBead::deprotonate()
    {
        protonationState_ = 0;
    }
    
    bool ContinuousProtonatableBead::isProtonated() const
    {
        return this->protonationState() == 1;
    }
    
    std::size_t ContinuousProtonatableBead::protonationState() const 
    {
        return protonationState_;
    }
    
    void ContinuousProtonatableBead::write(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        stream.setf(std::ios::scientific);        
        Particle::write(stream);
        stream << space << conf::CONTINUOUSLY_PROTONATABLE << space << protonationState_ 
               << space << std::setw(conf::WIDTH) << x_ 
               << space << std::setw(conf::WIDTH) << I_;
    }
    
    void ContinuousProtonatableBead::writeState(std::ostream& stream) const
    {
        const auto space = conf::SPACE;
        
        Particle::writeState(stream);
        stream << space << protonationState_ << space << x_ << space << I_;
    }
    
    void ContinuousProtonatableBead::readState(std::istream& stream)
    {
        Particle::readState(stream);
        stream >> protonationState_ >> x_ >> I_;        
    }
        
    cprot_bead_ptr_t ContinuousProtonatableBead::create(std::size_t id,
                                                        std::size_t index, 
                                                        const std::string& name,
                                                        std::size_t protonationState,
                                                        const spec_ptr_t& spec)
    {  
        return cprot_bead_ptr_t(
                new ContinuousProtonatableBead(id, 
                                               index, 
                                               name, 
                                               protonationState, 
                                               spec)
                );
    }
    
    ContinuousProtonatableBead::ContinuousProtonatableBead(std::size_t id,
                                                           std::size_t index, 
                                                           const std::string &name,
                                                           std::size_t protonationState,
                                                           const spec_ptr_t &spec) :
        Bead(id, index, name, spec), 
        x_{0.0}, I_{0.0}, protonationState_{protonationState}
    {        
        if ( protonationState_ > 1 ) {
            throw std::domain_error("ProtonatingBead: Illegal protonation state.");
        }        
    }
        
    std::ostream& operator << (std::ostream& stream, 
                               const ContinuousProtonatableBead& bead)
    {
        bead.write(stream);
        return stream;
    }
    
}