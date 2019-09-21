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
 * File:   properties.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 4, 2019, 3:42 PM
 */

#ifndef PPROPERTIES_HPP
#define PPROPERTIES_HPP

#include "ptypes.hpp"
#include "atom.hpp"
#include "bead.hpp"
#include <vector>

namespace simploce {
    namespace properties {
        
        /**
         * Finds particle in a given particle collection.
         * @param id Particle identifier.
         * @param particles Particles.
         * @return Particle or nullptr if the particle cannot be identified.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        std::shared_ptr<P> find(std::size_t id, const CONT<std::shared_ptr<P>>& particles)
        {
            for (auto p : particles) {
                if ( p->id() == id ) {
                    return p;
                }
            }
            return nullptr;
        }                    
        
        /**
         * Returns total charge of a collection of particles.
         * @param P Particle type.
         * @param particles Particles.
         * @return Total charge.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        mass_t mass(const CONT<std::shared_ptr<P>>& particles)
        {
            mass_t total{0.0};
            for (auto p : particles) {
                total += p->mass();
            }
            return total;
        }
           
        /**
         * Returns total charge of a collection of particles.
         * @param P Particle type.
         * @param particles Particles.
         * @return Total charge.
         */
        template <typename P, template<typename, typename ...> class CONT = std::vector>
        charge_t charge(const CONT<std::shared_ptr<P>>& particles)
        {
            charge_t total{0.0};
            for (auto p : particles) {
                total += p->charge();
            }
            return total;
        }
        
        /**
         * Returns the center of mass of a collection of particles.
         * @param P Particle type.
         * @param particles Particles.
         * @return Center of mass.
         */
        template <typename P, template <typename, typename...> class CONT = std::vector> 
        position_t centerOfMass(const CONT<std::shared_ptr<P>>& particles)
        {
            mass_t total{0.0};
            position_t cm{};
            for (auto p : particles) {
                total += p->mass();
                cm += p->mass()() * p->position();
            }
            return position_t{cm / total};
        }
    }
}

#endif /* PROPERTIES_HPP */

