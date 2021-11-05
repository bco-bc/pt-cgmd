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
 * File:   mc.hpp
 * Author: juffer
 *
 * Created on 19 October 2019, 18:56
 */

#ifndef MC_HPP
#define MC_HPP

#include "s-types.hpp"
#include "simploce/particle/bead.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Monte Carlo
     */
    template <typename P>
    class MC;
    
    /**
     * Specialization for beads.
     */
    template <>
    class MC<Bead> {
    public:
                
        /**
         * Constructor.
         * @param sm Simulation model. 
         */
        MC(const cg_sim_model_ptr_t& sm);
        
        /**
         * Performs the simulation.
         * @param param Parameters. Must provide,
         * <ul>
         *  <li>nsteps: Number of steps.</li>
         *  <li>
         *      nwrite: Number of steps between writing simulation data and saving state
         *      in the trajectectory.
         *  </li>
         * </ul>
         * @param trajStream Output trajectory stream.
         * @param dataStream Output simulation data stream.
         */
        void perform(const sim_param_t& param,
                     std::ofstream& trajStream,
                     std::ofstream& dataStream);
        
    private:
    
        cg_sim_model_ptr_t sm_;
        
    };
}

#endif /* MC_HPP */

