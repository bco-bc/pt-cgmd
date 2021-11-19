/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:36 PM
 */

#ifndef S_UTIL_HPP
#define S_UTIL_HPP

#include "s-types.hpp"
#include "simploce/units/units-mu.hpp"
#include "pair-lists.hpp"
#include "s-conf.hpp"
#include "simploce/particle/p-properties.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <vector>
#include <set>
#include <thread>

namespace simploce {
    namespace properties {

        /**
         * Returns inverse Debye screening length, in nm^-1.
         * @tparam T Particle type.
         * @param particles Particles.
         * @return Kappa.
         */
        template <typename T>
        real_t kappa(const std::vector<std::shared_ptr<T>>& particles) {
            return 0.0;
        }
        
        /**
         * Calculates instantaneous temperature for a collection of particles.
         * @param particles Particles.
         * @param eKin Kinetic energy.
         * @return Instantaneous temperature.
         */
        template <typename T>
        temperature_t temperature(const std::vector<std::shared_ptr<T>>& particles, 
                                  const energy_t& eKin)
        {
            std::size_t nParticles = particles.size();
            real_t nDof = 3.0 * real_t(nParticles) - 3.0;  // Assuming total momentum is constant.
            if (nDof > 3 ) {
                return 2.0 * eKin() / (nDof * units::mu<real_t>::KB );  // In K.
            } else {
                // No point calculating temperature for a low number of degrees of freedom.
                return 0.0;
            }        
        }
        
        /**
         * Returns pressure. Calculated from Virial Theorem.
         * @param particles All particles.
         * @param temperature Temperature.
         * @param box Simulation box.
         * @return Pressure.
         */
        template <typename T>
        pressure_t pressure(const std::vector<std::shared_ptr<T>>& particles,
                            const temperature_t& temperature,
                            const box_ptr_t& box)
        {
            volume_t volume = box->volume();
            real_t virial1 = 0.0;
            pressure_t pressure{};

            std::size_t nParticles = particles.size();
            for (const auto& particle : particles) {
                position_t r = particle->position();
                force_t f = particle->force();
                virial1 += inner<real_t>(f,r);
            }
            if ( volume() > 0.0 ) {
                virial1 /= ( 3.0 * volume() );
                real_t virial2 =
                        real_t(nParticles) * units::mu<real_t>::KB * temperature() / volume();
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
            static util::Logger logger("simploce::properties::tooClose");
            static length_t rMin = conf::SHORT;
        
            auto Rij = std::get<2>(ef);
            if (Rij() < rMin() ) {
                std::string message =
                        "WARNING: Rij < " + util::toString(rMin()) +
                        ", Rij = " + util::toString(Rij) +
                        ", pi = " + pi->name() + ", index = " + util::toString(pi->index()) +
                        ", id = " + util::toString(pi->id()) +
                        ", pj = " + pj->name() << ", index = " + util::toString(pj->index()) +
                        ", id = " + util::toString(pj->id()) +
                        ", energy: " + util::toString(std::get<0>(ef));
                logger.warn(message);
            }
        }
        
    }
}

#endif /* S_UTIL_HPP */

