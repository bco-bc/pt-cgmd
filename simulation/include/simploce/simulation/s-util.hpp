/*
 * Author: André H. Juffer.
 * Created on 15/06/2022.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_S_UTIL_HPP
#define SIMULATION_S_UTIL_HPP

#include "simploce/particle/p-util.hpp"

namespace simploce {
    namespace util {

        /**
         * Scales or adjust all velocities to a reference temperature. Use only for equilibration purposes.
         * @param particleSystem Particle system.
         * @param referenceTemperature Reference referenceTemperature.
         * @param isMesoscale If true, referenceTemperature is expressed as kT=n, where k is the
         * Boltzmann constant and n is a nonnegative integer.
         */
        void scaleVelocities(const p_system_ptr_t& particleSystem,
                             const temperature_t& referenceTemperature,
                             bool isMesoscale = false);

        /**
         * Removes overall center of mass motion (sets to zero).
         * @param particleSystem Particle system.
         */
        void removeCenterOfMassMotion(const p_system_ptr_t& particleSystem);
    }
}

#endif //SIMULATION_S_UTIL_HPP
