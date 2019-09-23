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
 * File:   analyzer.hpp
 * Author: juffer
 *
 * Created on 22 September 2019, 15:32
 */

#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include "atypes.hpp"

namespace simploce {
    
    /**
     * Performs an analysis on a collection of particles.
     * @param P Particle type.
     */
    template <typename P>
    struct Analyzer {
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        ~Analyzer() {}
        
        /**
         * Carries out the analysis.
         * @param particles Particles.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        virtual void perform(const std::vector<p_ptr_t> &particles) = 0;
    };
}

#endif /* ANALYZER_HPP */

