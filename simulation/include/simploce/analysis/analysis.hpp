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
 * File:   analasis.hpp
 * Author: juffer
 *
 * Created on 22 September 2019, 15:39
 */

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "atypes.hpp"
#include "simploce/particle/bead.hpp"
#include <iostream>

namespace simploce {

    /**
     * Specifies an analysis of a given simulation model.
     */
    template <typename P>
    class Analysis;
    
    template <>
    class Analysis<Bead> {
    public:
        
        /**
         * Constructor
         * @param sm Simulation model
         * @param analyzer Analyzer.
         */
        Analysis(const cg_sim_model_ptr_t& sm,
                 const cg_analyzer_ptr_t& analyzer);
        
        /**
         * Drives the analysis.
         * @param trajectory Input trajectory stream.
         */
        void perform(std::istream& trajectory);
        
    private:
        
        cg_sim_model_ptr_t sm_;
        cg_analyzer_ptr_t analyzer_;
    };    
    
}

#endif /* ANALYSIS_HPP */

