/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mu-units.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 3:40 PM
 */

#ifndef MU_UNITS_HPP
#define MU_UNITS_HPP

#include "si-units.hpp"

namespace simploce {
    
    /**
     * Numerical values of physical constants in 'molecular units' (MU units). 
     * These are derived from values in SI units. 
     * Units for MU are:
     * <ul>
     *   <li>Time: ps</li>
     *   <li>Distance, position: nm (= 10^-9 m)</li>
     *   <li>Velocity: nm/ps (1 ps = 10^-12 s)</li>
     *   <li>Mass: u (unified atomic mass unit, 1 u = 1.66054e-27 kg)</li>
     *   <li>Momentum: (u nm)/ps</li>
     *   <li>Energy: kJ/mol = (u nm^2)/(ps^2)</li>
     *   <li>Force: kJ/(mol nm) = (u nm)/(ps^2)</li>
     *   <li>Charge: e (= 1.6021766208e-19 C).</li>
     * </ul>
    */
    template <typename V>
    struct MUUnits {
    
        /**
         * Elementary charge, in e.
         */
        static const V E;      
    
        /**
         * Electric constant. In (mol e^2)/(kJ nm)
         */
        static const V E0;
    
        /**
         * 4 * PI * E0. In (kJ nm) / (mol e2).
         */
        static const V FOUR_PI_E0;
    
        /**
         * Electric factor 1.0/(4 * PI * E0). In (kJ nm)/(mol e^2)
         */
        static const V F_EL;
    
        /**
         * Boltzmann constant. In kJ/(mol K).
         */
        static const V KB;

        /**
         * Molar gas constant. In kJ/(mol K).
         */
        static const V R;
    
        /**
         * Faraday constant. In e/mol.
         */    
        static const V F;
    
        /**
         * Value of kT at room temperature. In kJ/mol.
         */
        static const V kT;
    
        /**
         * Proton mass. In u.
         */
        static const V PROTON_MASS;

        /**
         * Proton charge. In e.
         */
        static const V PROTON_CHARGE;

        /**
         * Viscosity of water at 298.15 K. In u/(nm ps).
         */
        static const V WATER_VISCOSITY;
    
        /**
         * l to nm^3 conversion factor.
         */
        static const V l_to_nm3;

        /**
         * Ångström to nm conversion factor.
         */
        static const V Angstrom_to_nm;

        /**
         * V to kJ/(mol e).
         */
        static const V V_to_kJ_mol_e;
        
    };
  
    template <typename V>
    const V MUUnits<V>::E = 1.0;
  
    template <typename V>
    const V MUUnits<V>::E0 = 
        SIUnits<V>::E0/(SIUnits<V>::E * 1.0e+09) * 1.0e+03 / (SIUnits<V>::eV * SIUnits<V>::NA);
  
    template <typename V>
    const V MUUnits<V>::FOUR_PI_E0 = 4.0 * MathConstants<V>::PI * MUUnits<V>::E0;
  
    template <typename V>
    const V MUUnits<V>::F_EL = 1.0 / MUUnits<V>::FOUR_PI_E0;
  
    template <typename V>
    const V MUUnits<V>::KB = SIUnits<V>::KB * SIUnits<V>::NA / 1.0e+03;

    template <typename V>
    const V MUUnits<V>::R = SIUnits<V>::R / 1.0e+03;

    template <typename V>
    const V MUUnits<V>::F = SIUnits<V>::F / SIUnits<V>::E;
  
    template <typename V>
    const V MUUnits<V>::kT = MUUnits<V>::KB * SIUnits<V>::roomT;
  
    template <typename V>
    const V MUUnits<V>::PROTON_MASS = SIUnits<V>::PROTON_MASS / SIUnits<V>::MU;

    template <typename V>
    const V MUUnits<V>::PROTON_CHARGE = 1.0;

    template <typename V>
    const V MUUnits<V>::WATER_VISCOSITY =
        SIUnits<V>::WATER_VISCOSITY / (SIUnits<V>::MU * 1.0e+09 * 1.0e+12);
  
    template <typename V>
    const V MUUnits<V>::l_to_nm3 = 1.0e-03 * 1.0e+27;

    template <typename V>
    const V MUUnits<V>::Angstrom_to_nm = 0.1;

    template <typename V>
    const V MUUnits<V>::V_to_kJ_mol_e = SIUnits<V>::NA * SIUnits<V>::E/1.0e+03;
      
}

#endif /* MU_UNITS_HPP */

