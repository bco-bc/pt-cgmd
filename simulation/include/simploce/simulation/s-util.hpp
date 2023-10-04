/*
 * Author: André H. Juffer.
 * Created on 15/06/2022.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_S_UTIL_HPP
#define SIMULATION_S_UTIL_HPP

#include "simploce/particle/p-util.hpp"
#include "simploce/types/s-types.hpp"

namespace simploce {
    namespace util {

        /**
         * Returns cutoff distance for pair list calculation.
         * @param param Simulation parameters. The cutoff for long-range interactions is assumed here.
         * @param particleSystem Particle system.
         * @return Cutoff distance for pair lists.
         */
        dist_t computePairListCutoff(const param_ptr_t& param,
                                      const p_system_ptr_t & particleSystem);

        /**
         * Scales or adjust all velocities to a reference temperature. Use only for equilibration purposes.
         * @param all Particle system.
         * @param referenceTemperature Reference referenceTemperature.
         * @param mesoscopic If true, referenceTemperature is expressed as kT=n, where k is the
         * Boltzmann constant and n is a nonnegative integer.
         */
        void scaleVelocities(const p_system_ptr_t& all,
                             const temperature_t& referenceTemperature,
                             bool mesoscopic = false);

        /**
         * Removes overall center of mass motion (sets to zero).
         * @param particleSystem Particle system.
         */
        void removeCenterOfMassMotion(const p_system_ptr_t& particleSystem);

        /**
         * Place particles inside the simulation box.
         * @param particleSystem Particle system.
         * @param bc Boundary conditions.
         */
        void placeInsideBox(const p_system_ptr_t& particleSystem, const bc_ptr_t& bc);
    }
}

#endif //SIMULATION_S_UTIL_HPP
