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
#include "sconf.hpp"
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
         * Returns pressure.
         * @param all All particles.
         * @param temperature Temperature.
         * @param box Simulation box.
         * @return Pressure.
         */
        template <typename T>
        pressure_t pressure(const std::vector<std::shared_ptr<T>>& all,
                            const temperature_t& temperature,
                            const box_ptr_t& box)
        {
            volume_t volume = box->volume();
            mass_t totalMass = 0.0;
            real_t virial1 = 0.0;
            position_t centerOfMass{};
            momentum_t totalMomentum{};
            pressure_t pressure{};

            std::size_t nparticles = all.size();
            for (const auto& particle : all) {
                mass_t mass = particle->mass()();
                position_t r = particle->position();
                momentum_t p = particle->momentum();
                force_t f = particle->force();
                totalMass += mass;
                centerOfMass += ( mass() * r );
                totalMomentum += p;
                virial1 += inner<real_t>(f,r);
            }
            centerOfMass /= totalMass();
            if ( volume() > 0.0 ) {
                virial1 /= ( 3.0 * volume() );
                real_t virial2 = nparticles * MUUnits<real_t>::KB * temperature() / volume();
                pressure = virial1 + virial2; // In kJ/(mol nm^3)    
            } else {
                pressure = 0.0;
            }
            return pressure;            
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
         * Returns dielectric constant according to Fröhlich.
         * @param aveM2 The average of the M*M, where M is the total dipole moment.
         * @param temperature Temperature.
         * @param box Simulation box.
         * @return Number.
         * @see 
         */
        real_t frohlich(real_t aveM2, 
                        const temperature_t& temperature,
                        const box_ptr_t& box);
        
        /**
         * Writes a warning to std::clog if particles are too close.
         * @param pi Particle 1
         * @param pj Particle 2
         * @param ef Holds energy and forces.
         */
        template <typename P>
        static void 
        tooClose(const std::shared_ptr<P>& pi, 
                 const std::shared_ptr<P>& pj, 
                 const std::tuple<energy_t, force_t, length_t>& ef)
        {
            static length_t Rmin = conf::CLOSE;
        
            auto Rij = std::get<2>(ef);
            if ( Rij() < Rmin() ) {
                std::clog << "WARNING: Rij < " << Rmin() << ", Rij = " << Rij 
                          << " pi = " << pi->name() << ", index = " << pi->index()
                          << " pj = " << pj->name() << ", index = " << pj->index()
                        << " energy: " << std::get<0>(ef) 
                        << std::endl;            
            }
        }
        
    }
}

#endif /* SIM_UTIL_HPP */

