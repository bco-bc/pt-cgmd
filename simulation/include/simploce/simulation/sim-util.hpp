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
 * File:   sim-util.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:36 PM
 */

#ifndef SIM_UTIL_HPP
#define SIM_UTIL_HPP

#include "stypes.hpp"
#include "simploce/util/mu-units.hpp"
#include "pair-lists.hpp"
#include <vector>
#include <set>
#include <thread>

namespace simploce {
    namespace util {        
        
        /**
         * Calculates instantaneous temperature for a collection of particles.
         * @param particles Particles.
         * @param ekin Kinetic energy.
         * @return Instantaneous temperature.
         */
        template <typename T>
        temperature_t temperature(const std::vector<std::shared_ptr<T>>& particles, 
                                  const energy_t& ekin)
        {
            std::size_t nparticles = particles.size();
            //std::size_t ndof = ( 3 * nparticles - 3 );  // -3 to remove rigid body translation
            std::size_t ndof = 3 * nparticles;
            if ( ndof > 3 ) {
                return 2.0 * ekin() / ( real_t(ndof) * MUUnits<real_t>::KB );  // In K.
            } else {
                // No point calculating temperature for a low number of degrees of freedom.
                return 0.0;
            }        
        }            
        
        /**
         * Returns cutoff distance.
         * @param box Simulation box.
         * @return Cutoff distance. Always <= 0.5 * box.size().
         */
        length_t cutoffDistance(const box_ptr_t& box);
        
        /**
         * Returns square of cutoff distance.
         * @param box Simulation box.
         * @return Square of cutoff distance.
         */
        real_t squareCutoffDistance(const box_ptr_t& box);
        
        /**
         * Returns sub lists of a list of items.
         * @param items Items
         * @return Sub lists of items.
         */
        template <typename T>
        std::set<std::vector<T>> 
        makeSubLists(const std::vector<T>& items)
        {
            using sublists_t = std::set<std::vector<T>>;
            
            if ( items.empty() ) {
                return sublists_t{};
            }
            
            // Sub lists.
            sublists_t subLists;
            
            std::size_t counter = 0;      
            static const std::size_t nsublists = std::thread::hardware_concurrency();        
            std::size_t numberOfItemsPerSubList = items.size() / nsublists;                                                              
            for (std::size_t k = 0; k != nsublists; ++k) {
                std::vector<T> single{};  // One sublist of items.
                std::size_t n = 0;
                while (counter != items.size() && n != numberOfItemsPerSubList) {
                    single.push_back(items[counter]);
                    counter += 1;
                    n += 1;
                }
                subLists.push_back(single);          // Store list.
            }
    
            // Add remaining pairs, if any, to the last set.
            if ( counter < items.size() ) {
                auto& last = *(items.end() - 1);
                while ( counter != items.size() ) {
                    last.push_back(items[counter]);
                    counter += 1;
                }
            }
            
            return subLists;
        }
        
    }
}

#endif /* SIM_UTIL_HPP */

