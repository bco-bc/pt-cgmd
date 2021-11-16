/*
 * File:   sim-util.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:38 PM
 */

#include "simploce/simulation/s-properties.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/units/units-mu.hpp"

namespace simploce {
    namespace properties {
        
        length_t cutoffDistance(const box_ptr_t& box)
        {
            length_t rc = 0.5 * box->size();
            return conf::CUTOFF_DISTANCE() > rc() ? rc : conf::CUTOFF_DISTANCE;
        }
        
        real_t squareCutoffDistance(const box_ptr_t& box)
        {
            length_t rc = cutoffDistance(box);
            return rc() * rc();
        }
        
        real_t frohlich(real_t aveM2, 
                        const temperature_t& temperature,
                        const box_ptr_t& box)
        {
            auto E0 = units::mu<real_t>::E0;
            auto kT = units::mu<real_t>::KB * temperature();
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