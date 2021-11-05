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
 * File:   constant-rate-pt.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 4, 2019, 2:16 PM
 */

#include "simploce/simulation/pt-pair-list-generator.hpp"
#include "simploce/simulation/continuous.hpp"
#include "simploce/simulation/constant-rate-pt.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/particle/protonatable-bead.hpp"
#include "simploce/util/telegraph-process.hpp"
#include <stdexcept>

namespace simploce {
    
    static const rate_t RATE = 1.0 / 1.5;
    static const real_t GAMMA = 1.0 / RATE(); 
    
    using prot_pair_t = ProtonTransferPairListGenerator::prot_pair_t;
    
    static bool hasProton_(const prot_pair_t& pair)
    {
        return pair.first->isProtonated() || pair.second->isProtonated();
    }
    
    ConstantRateProtonTransfer::ConstantRateProtonTransfer() :
        rate_(RATE), gamma_(GAMMA)
    {        
    }
    
    ConstantRateProtonTransfer::ConstantRateProtonTransfer(const rate_t& rate, 
                                                           const real_t& gamma) :
        rate_{rate}, gamma_{gamma}
    {   
        if ( rate_() < 0.0 ) {
            throw std::domain_error(
                "ConstantRateProtonTransfer: "
                "Rate value must be a nonnegative number."
            );
        }
        if ( gamma_ < 0.0 ) {
            throw std::domain_error(
                "ConstantRateProtonTransfer: "
                "Inverse of decay time must be a nonnegative number."
            );
        }
    }   
        
    void 
    ConstantRateProtonTransfer::transfer(const sim_param_t& param,
                                         const std::vector<prot_bead_ptr_t>& continuous,
                                         const prot_pair_list_t& pairList) const
    {
        stime_t dt = param.get<real_t>("timestep");
        
        // Current protonation states.
        std::vector<std::size_t> dX(continuous.size(), 0);

        // Protonation states.
        velocity_t u{0, 0, 0};
        std::size_t npairs = 0;
        for (auto pair : pairList) {
            if ( hasProton_(pair) ) {
                npairs += 1;
                auto p1 = pair.first;
                std::size_t index1 = p1->index();
                auto state_1 = p1->protonationState();
                int X1 = state_1.isProtonated() ? 1 : 0;
                auto p2 = pair.second;
                std::size_t index2 = p2->index();
                auto state_2 = p2->protonationState();
                int X2 = state_2.isProtonated() ? 1 : 0;
                
                // Transfer only if one of the beads is protonated.
                int dX1 = 0;
                int dX2 = 0;
                if ( X1 == 1 && X2 == 0) {
            
                    // From X to Y.
                    //int T = 2 * X - 1;
                    int dT = TelegraphProcess::increment(X1, dt, rate_);
                    dX1 = dT / 2;
                    dX2 = -dX1;
                    X1 += dX1;
                    X2 += dX2;
                    if ( X1 == 0 ) {
                        state_1.protonate();
                    }
                    if ( X2 == 1) {
                        state_2.deprotonate();
                    }
                } else {
                    if ( X1 == 0 && X2 == 1 ) {
              
                        // From Y to X.
                        //int T = 2 * Y - 1;
                        int dT = TelegraphProcess::increment(X2, dt, rate_);
                        dX2 = dT / 2;
                        dX1 = -dX2;
                        X1 += dX1;
                        X2 += dX2;
                        if ( X1 == 1 ) {
                            state_1.protonate();
                        }
                        if ( X2 == 0 ) {
                            state_2.deprotonate();
                        }
                    }
                }
                
                dX[index1] += dX1;
                dX[index2] += dX2;
                
                u += (p1->velocity() + p2->velocity());
            }
        }
        
        // Assume an average.
        u *= 0.5 / real_t(npairs);
        
        for (auto p : continuous) {
            
            mass_t m = p->mass();  // Just before mass transfer.
            
            /* Update state x and current I of protonatable.
            std::size_t index = p->index();
            real_t x = p->state();
            real_t I = p->current();
            real_t dI = -gamma_ * I + gamma_ * dX[index];
            I += dI;
            real_t dx = I * dt();
            x += dx;
            p->state(x);
            p->current(I);
            
            // Update momentum due to mass transfer.
            mass_t dm = dI * conf::MASS_PROTON;
            velocity_t v = p->velocity();
            static velocity_t dv;
            for (std::size_t k = 0; k != 3; ++k) {
                dv[k] = - (v[k] - u[k]) * dm() / m();
            }
            v += dv;
            p->velocity(v);
             */
        }
        
        // Done.
    }
    
}