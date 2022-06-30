/*
 * Created by André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 11/29/21.
 */

#ifndef SIMULATION_SP_CONF_HPP
#define SIMULATION_SP_CONF_HPP

#include <string>

namespace simploce {
    namespace conf {

        // Pair potential types.
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
        const std::string LJ_SF = "lj+sf";      // Lennard-Jones plus shifted force electrostatics.
        const std::string HP = "hp";            // Harmonic.
        const std::string HA_QP = "ha-qp";      // Halve-attractive quartic.
        const std::string SR = "sr";            // Soft repulsive (for MVV_DPD).
        const std::string HP_SR = "hp+sr";      // Harmonic + soft repulsive.

        // External potential types.
        const std::string C_SRF_CG_DENSITY = "c-srf-cg-dens";            // Const surface charge density.
        const std::string ELECTRIC_POTENTIAL_DIFFERENCE = "epd";         // Electric potential difference.
        const std::string WALL = "wall";                                 // Impenetrable wall.
    }
}

#endif //SIMULATION_SP_CONF_HPP
