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
        
        const std::string LEAP_FROG = "lf";
        const std::string LANGEVIN_VELOCITY_VERLET = "lvv";
        const std::string PT_LANGEVIN_VELOCITY_VERLET = "pt-lvv";
        const std::string VELOCITY_VERLET = "vv";
        
        const std::string NOBC = "no-bc";
        const std::string PBC = "pbc";
        
        const std::string POLARIZABLE_WATER = "pol-water";
        const std::string ACID_BASE_SOLUTION = "acid-base-solution";
        const std::string ELECTROLYTE = "electrolyte";
        const std::string LJ_FLUID = "lj-fluid";
        const std::string HP = "hp";
        
        // Default cutoff distance for non-bonded interactions.
        static length_t CUTOFF_DISTANCE{2.6};  // nm.
        
        // Minimum number of particles.
        const std::size_t MIN_NUMBER_OF_PARTICLES = 1000;    
        
        /**
         * Default cutoff distance beyond which protonatables cannot transfer protons.
         */
        static length_t CUTOFF_DISTANCE_PT{0.4};

        /**
         * Minimum value below which a distance is considered to be close.
         */
        static length_t CLOSE{0.15};
    }
}

#endif /* SCONF_HPP */

