/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:36 PM
 */

#ifndef SIM_UTIL_HPP
#define SIM_UTIL_HPP

#include "s-types.hpp"
#include "simploce/units/units-mu.hpp"
#include "pair-lists.hpp"
#include "s-conf.hpp"
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
            real_t ndof = 3 * nparticles - 3;  // Assuming total momentum is constant.
            if ( ndof > 3 ) {
                return 2.0 * ekin() / ( ndof * units::mu<real_t>::KB );  // In K.
            } else {
                // No point calculating temperature for a low number of degrees of freedom.
                return 0.0;
            }        
        }
        
        /**
         * Returns pressure. Calculated from Virial Theorem.
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
            real_t virial1 = 0.0;
            pressure_t pressure{};

            std::size_t nparticles = all.size();
            for (const auto& particle : all) {
                position_t r = particle->position();
                force_t f = particle->force();
                virial1 += inner<real_t>(f,r);
            }
            if ( volume() > 0.0 ) {
                virial1 /= ( 3.0 * volume() );
                real_t virial2 = 
                    nparticles * units::mu<real_t>::KB * temperature() / volume();
                pressure = virial2 - virial1; // In kJ/(mol nm^3)    
            } else {
                pressure = 0.0;
            }
            return pressure;            
        }
        
        /**
         * Returns cutoff distance.
         * @param box Simulation box.
         * @return Cutoff distance. Always <= 0.5 * box.size().
         * @see conf::CUTOFF_DISTANCE
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
                          << ", pi = " << pi->name() << ", index = " << pi->index()
                          << ", id = " << pi->id() 
                          << ", pj = " << pj->name() << ", index = " << pj->index()
                          << ", id = " << pj->id()
                          << ", energy: " << std::get<0>(ef) 
                          << std::endl;            
            }
        }
        
    }
}

#endif /* SIM_UTIL_HPP */

