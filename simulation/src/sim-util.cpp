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
 * File:   sim-util.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:38 PM
 */

#include "simploce/simulation/sim-util.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/util/mu-units.hpp"

namespace simploce {
    namespace util {
        
        length_t cutoffDistance(const box_ptr_t& box)
        {
            return 0.5 * box->size();            
            /*
            length_t halve = 0.5 * box->size();
            return conf::RCUTOFF_DISTANCE_() > halve() ? 
                   halve : 
                   conf::RCUTOFF_DISTANCE_;             
             */
        }
        
        real_t squareCutoffDistance(const box_ptr_t& box)
        {
            length_t rc = cutoffDistance(box);
            return rc * rc;
        }
        
        real_t frohlich(real_t aveM2, 
                        const temperature_t& temperature,
                        const box_ptr_t& box)
        {
            auto E0 = MUUnits<real_t>::E0;
            auto kT = MUUnits<real_t>::KB * temperature();
            auto volume = box->volume();
            real_t h = 1.0 / (E0 * volume) * aveM2 / (3.0 * kT);
            real_t a = 2;
            real_t b = -1 - 3 * h;
            real_t c = -1;
            real_t D = b * b - 4 * a * c;
            assert(D > 0);
            return (-b + std::sqrt(D)) / (2.0 * a);
        }
                
    }
}