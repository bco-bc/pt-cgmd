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
        const std::string PBC = "pbc";       // Periodic boundary conditions.

        // Algorithm.
        const std::string LEAP_FROG = "lf";
        const std::string LANGEVIN_VELOCITY_VERLET = "lvv";
        const std::string PT_LANGEVIN_VELOCITY_VERLET = "pt-lvv";    // With proton transfer (PT).
        const std::string VELOCITY_VERLET = "vv";
        
        const std::string POLARIZABLE_WATER = "pol-water";
        const std::string ACID_BASE_SOLUTION = "acid-base-solution";
        const std::string ELECTROLYTE = "electrolyte";
        const std::string LJ_FLUID = "lj-fluid";

        // Interaction type.
        const std::string LJ = "lj";             // Lennard-Jones
        const std::string LJ_RF = "lj+rf";       // Lennard-Jones + reaction field.
        const std::string COULOMB = "coulomb";   // Coulomb.
        const std::string HP = "hp";             // Harmonic.
        const std::string HA_QP = "ha-qp";       // Halve-attractive quartic.
        
        // Default cutoff distance for non-bonded interactions.
        const length_t CUTOFF_DISTANCE{2.6};  // nm.
        
        // Minimum number of particles.
        const std::size_t MIN_NUMBER_OF_PARTICLES{500};
        
        /**
         * Default cutoff distance beyond which protonatables cannot transfer protons.
         */
        static length_t CUTOFF_DISTANCE_PT{0.4};

        /**
         * Value below which a distance is considered to be close.
         */
        static distance_t CLOSE{0.15};
    }
}

#endif /* SCONF_HPP */

