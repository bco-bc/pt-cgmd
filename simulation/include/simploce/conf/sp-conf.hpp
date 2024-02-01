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

        // Non-bonded/bonded pair potentials.
        const std::string LJ = "lj";                                     // Lennard-Jones
        const std::string LJ_RF = "lj+rf";                               // Lennard-Jones + reaction field electrostatics.
        const std::string RF = "rf";                                     // Reaction field electrostatics.
        const std::string HS_SF = "hs+sf";                               // Hard sphere potential +
                                                                         // shifted force electrostatics
                                                                         // NOTE: Must only be used in combination with
                                                                         // Monte Carlo.
        const std::string HS_SC = "hs+sc";                               // Hard sphere potential + screened Coulomb
                                                                         // electrostatics
                                                                         // NOTE: Must only be used in combination
                                                                         // with Monte Carlo.
        const std::string HS_RF = "hs+rf";                               // Hard sphere potential + reaction field
                                                                         // electrostatics.
                                                                         // Must only be used in combination with
                                                                         // Monte Carlo.
        const std::string SC = "sc";                                     // Screened Coulomb electrostatics.
        const std::string SF = "sf";                                     // Shifted force electrostatics.
        const std::string LJ_SF = "lj+sf";                               // Lennard-Jones plus shifted force
                                                                         // electrostatics.
        const std::string SR = "sr";                                     // Soft repulsive.
        const std::string HP = "hp";                                     // Harmonic.
        const std::string HA_QP = "ha-qp";                               // Halve-attractive quartic.
        const std::string HA_HP = "ha-hp";                               // Halve-attractive harmonic.
        const std::string HP_SR = "hp+sr";                               // Harmonic + soft repulsive.
        const std::string GAUSS_SF = "g_sf";                             // Overlapping Gaussian charge densities with
                                                                         // shift force electrostatics.
        const std::string GAUSS_SF_SR = "g_sf+sr";                       // Overlapping Gaussian charge densities with
                                                                         // shift force electrostatics +
                                                                         // soft repulsive.

        const std::string NONE_INTERACTING = "none";                     // None-interacting potential.

        // External potential types.
        // const std::string C_SRF_CG_DENSITY = "c-srf-cg-dens";         // Constant surface charge density.
        const std::string VOLTAGE = "voltage";                           // Static uniform external electric field
                                                                         // in z-direction
        const std::string WALL = "wall";                                 // Impenetrable flat wall.
        const std::string PRESSURE_GRADIENT = "pres_grad";               // Pressure gradient in z-direction.
    }
}

#endif //SIMULATION_SP_CONF_HPP
