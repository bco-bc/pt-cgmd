/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 3, 2019, 11:37 AM
 */

#ifndef SCONF_HPP
#define SCONF_HPP

#include "simploce/conf/p-conf.hpp"
#include "sp-conf.hpp"

namespace simploce {
    namespace conf {

        // Boundary conditions
        const std::string NO_BC = "no-bc";   // None, no boundary conditions.
        const std::string PBC = "pbc";       // Periodic boundary conditions.

        // Algorithm ("Displacers").
        const std::string LEAP_FROG = "lf";
        const std::string LANGEVIN_VELOCITY_VERLET = "lvv";
        const std::string PT_LANGEVIN_VELOCITY_VERLET = "pt-lvv";    // With proton transfer (PT).
        const std::string VELOCITY_VERLET = "vv";
        const std::string MONTE_CARLO = "mc";
        const std::string DPD = "mvv-dpd";
        const std::string S1_DPD = "s1-dpd";

        // Minimum number of particles for concurrent simulations.
        const std::size_t MIN_NUMBER_OF_PARTICLES{500};
        
        /**
         * Default cutoff distance beyond which protonatables cannot transfer protons, in nm.
         */
        const length_t PT_CUTOFF_DISTANCE{0.4};

        /**
         * Value below which a distance is considered to be "short", in nm.
         */
        const dist_t SHORT{0.15};

        /**
         * A very large larger positive number.
         */
        const real_t LARGE = 1.0e+30;

        /**
         * A very small positive real number.
         */
        const real_t SMALL = 1.0e-20;

   }
}

#endif /* SCONF_HPP */

