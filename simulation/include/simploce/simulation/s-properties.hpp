/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 1:36 PM
 */

#ifndef S_UTIL_HPP
#define S_UTIL_HPP

#include "simploce/types/s-types.hpp"
#include "simploce/particle/p-properties.hpp"
#include <vector>

namespace simploce {
    namespace properties {

        /**
         * Returns inverse Debye screening length, in nm^-1.
         * @tparam T Particle typeName.
         * @param particles Particles.
         * @return Kappa.
         */
        real_t kappa(const std::vector<p_ptr_t>& particles);
        
        /**
         * Calculates instantaneous (kinetic) temperature for a collection of particles. Only
         * displaceable (not frozen) particles contribute to the temperature.
         * @param particles Particles.
         * @param eKin Kinetic energy.
         * @param mesoscopic Particles represent a isMesoscale system. If true, than an unit system is assumed
         * identical to that of dissipative particle dynamics, such that kT=n where k is Boltzmann constant and T
         * is the temperature, and n is an non-negative number.
         * Otherwise molecular units are assumed.
         * @return Instantaneous temperature.
         */
        temperature_t kineticTemperature(const std::vector<p_ptr_t>& particles,
                                         const energy_t& eKin,
                                         bool mesoscopic = false);
        
        /**
         * Returns pressure. Calculated from Virial theorem. Only displaceable (not frozen)
         * particles contribute to the pressure.
         * @param particles All particles.
         * @param temperature Temperature.
         * @param box Simulation box.
         * @return Pressure.
         */
        pressure_t pressure(const std::vector<p_ptr_t>& particles,
                            const temperature_t& temperature,
                            const box_ptr_t& box,
                            bool mesoscopic = false);
        
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
         * Writes a warning if particles are too close.
         * @param pi Particle 1
         * @param pj Particle 2
         * @param efl Holds energy, forces, and distance between particles 1 and 2.
         */
        static void 
        tooClose(const p_ptr_t& pi,
                 const p_ptr_t& pj,
                 const std::tuple<energy_t, force_t, length_t>& efl);

    }
}

#endif /* S_UTIL_HPP */

