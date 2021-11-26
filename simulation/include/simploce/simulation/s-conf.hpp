/*
 * File:   sconf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 11:37 AM
 */

#ifndef SCONF_HPP
#define SCONF_HPP

#include "simploce/particle/p-conf.hpp"


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

        // Interaction types.
        const std::string LJ = "lj";            // Lennard-Jones
        const std::string LJ_RF = "lj+rf";      // Lennard-Jones + reaction field electrostatics.
        const std::string RF = "rf";            // Reaction field electrostatics.
        const std::string HS_SF = "hs+sf";      // Hard sphere potential + shifted force electrostatics
                                                // Must only be used in combination with Monte Carlo.
        const std::string HS_SC = "hs+sc";      // Hard sphere potential + screened Coulomb electrostatics
                                                // Must only be used in combination with Monte Carlo.
        const std::string HS_RF = "hs+rf";      // Hard sphere potential + reaction field electrostatics.
                                                // Must only be used in combination with Monte Carlo.
        const std::string SC = "sc";            // Screened Coulomb electrostatics.
        const std::string SF = "sf";            // Shifted force electrostatics.
        const std::string HP = "hp";            // Harmonic.
        const std::string HA_QP = "ha-qp";      // Halve-attractive quartic.

        // Default cutoff distance for non-bonded interactions, in nm.
        const distance_t CUTOFF_DISTANCE{2.6};
        
        // Minimum number of particles.
        const std::size_t MIN_NUMBER_OF_PARTICLES{500};
        
        /**
         * Default cutoff distance beyond which protonatables cannot transfer protons, in nm.
         */
        const length_t PT_CUTOFF_DISTANCE{0.4};

        /**
         * Value below which a distance is considered to be "short".
         */
        const distance_t SHORT{0.15};

        /**
         * A very large larger positive number.
         */
        const real_t LARGE = 1.0e+30;
    }
}

#endif /* SCONF_HPP */

