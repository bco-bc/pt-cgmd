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
 * File:   constant-rate-pt.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 4, 2019, 2:11 PM
 */

#ifndef CONSTANT_RATE_PT_HPP
#define CONSTANT_RATE_PT_HPP

#include "pt.hpp"

namespace simploce {
    
    /**
     * Transfer protons with the same constant proton transfer event rate.
     */
    class ConstantRateProtonTransfer : public ProtonTransfer {
    public:
        
        /**
         * Default constructor with default values for rate and gamma.
         */
        ConstantRateProtonTransfer();
        
        /**
         * Constructor
         * @param rate Rate.
         * @param gamma gamma^-1 is time constant of decay.
         */
        ConstantRateProtonTransfer(const rate_t& rate, 
                                   const real_t& gamma);
        
        void transfer(const sim_param_t& param,
                      const std::vector<prot_bead_ptr_t>& continuous,
                      const prot_pair_list_t& pairList) const override;
        
    private:
        
        rate_t rate_;
        real_t gamma_;
    };
}

#endif /* CONSTANT_RATE_PT_HPP */

