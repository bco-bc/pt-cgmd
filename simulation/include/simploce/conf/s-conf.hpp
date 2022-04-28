/*
 * File:   sconf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 11:37 AM
 */

#ifndef SCONF_HPP
#define SCONF_HPP

#include "simploce/particle/p-conf.hpp"
#include "simploce/potentials/sp-conf.hpp"

namespace simploce {
    namespace conf {

        // Boundary conditions
        const std::string NO_BC = "no-bc";   // None.
        const std::string PBC = "boundaryCondition";       // Periodic boundary conditions.

        // Algorithm ("Displacers").
        const std::string LEAP_FROG = "lf";
        const std::string LANGEVIN_VELOCITY_VERLET = "lvv";
        const std::string PT_LANGEVIN_VELOCITY_VERLET = "pt-lvv";    // With proton transfer (PT).
        const std::string VELOCITY_VERLET = "vv";
        const std::string MONTE_CARLO = "mc";

        // Minimum number of particles.
        const std::size_t MIN_NUMBER_OF_PARTICLES{500};
        
        /**
         * Default cutoff distance beyond which protonatables cannot transfer protons, in nm.
         */
        const length_t PT_CUTOFF_DISTANCE{0.4};

        /**
         * Value below which a distance is considered to be "short".
         */
        const dist_t SHORT{0.15};

        /**
         * A very large larger positive number.
         */
        const real_t LARGE = 1.0e+30;

        /**
         * A small positive real number.
         */
        const real_t SMALL = 1.0e-10;

        // Default values

   }
}

#endif /* SCONF_HPP */

