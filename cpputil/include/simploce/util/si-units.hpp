/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   si-units.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 3:34 PM
 */

#ifndef SI_UNITS_HPP
#define SI_UNITS_HPP

#include "math-constants.hpp"

namespace simploce {
    
    /**
     * Numerical values of physical constants, in SI units.
     * @see <a href="http://physics.nist.gov/cuu/Constants/index.html">NIST</a>
     * @param V value type (float, double).
     */    
    template <typename V>
    struct SIUnits {
    
        /**
         * Elementary charge, in C.
         */
        static const V E;
    
        /**
         * Electric constant or vacuum permittivity or permittivity of free space. In F/m = C/(V m).
         */
        static const V E0;
    
        /**
         * 4 * PI * E0, in C/(V m).
         */
        static const V FOUR_PI_E0;
    
        /**
         * Electric factor 1.0/(4 * PI * E0), in m/F = (V m)/C.
         */
        static const V F_EL;
    
        /**
         * Boltzmann constant. In J/K.
        */
        static const V KB;

        /**
         * Molar gas constant. In J/(mol K).
         */
        static const V R;

        /**
         * Faraday constant. In C/mol.
        */
        static const V F;
    
        /**
         * Room temperature (25 degrees Celsius, 298.15 K), in K.
         */
        static const V roomT;
    
        /**
         * Value of kT at room temperature. In J.
         */
        static const V kT;
    
        /**
         * Avogadro constant. In 1/mol.
         */
        const static V NA;
    
        /**
         * Unified atomic mass unit. In kg.
         */
        static const V MU;
        
        /**
         * Proton mass. In kg.
         */
        static const V PROTON_MASS;

        /**
         * Proton charge. In C.
         */
        static const V PROTON_CHARGE;
    
        /**
         * Fine structure constant. Dimensionless.
         */
        static const V ALPHA;
    
        /**
         * Planck constant. In Js.
         */
        static const V H;
    
        /**
         * Magnetic constant. In N/(A^2) = (V s)/(A m).
         */
        static const V MU_0;
    
        /**
         * Speed of light. In m/s.
         */
        static const V C_0;
    
        /**
         * Value of 1eV = sqrt(2*h*alpha / (mu_0 * c_0))/1C. In J.
         * Approximately 1.60218 x 10-19 J.
         * @see <a href="http://en.wikipedia.org/wiki/Electronvolt">Wikipedia</a>
        */
        static const V eV;

        /**
         * Viscosity of water at 298.15 K. In kg m^-1 s^−1.
        */
        static const V WATER_VISCOSITY;

    };
  
    template <typename V>
    const V SIUnits<V>::E = 1.6021766208e-19;
  
    template <typename V>
    const V SIUnits<V>::E0 = 8.854187817e-12;
  
    template <typename V>
    const V SIUnits<V>::FOUR_PI_E0 = 4.0 * MathConstants<V>::PI * SIUnits<V>::E0;
  
    template <typename V>
    const V SIUnits<V>::F_EL = 1.0 / SIUnits<V>::FOUR_PI_E0;
  
    template <typename V>
    const V SIUnits<V>::KB = 1.38064852e-23;

    template <typename V>
    const V SIUnits<V>::NA = 6.022140857e+23;

    template <typename V>
    const V SIUnits<V>::R = SIUnits<V>::KB * SIUnits<V>::NA;

    template <typename V>
    const V SIUnits<V>::F = SIUnits<V>::E * SIUnits<V>::NA;
  
    // see https://physics.nist.gov/cuu/Units/outside.html
    template <typename V>
    const V SIUnits<V>::MU = 1.66054e-27;
  
    template <typename V>
    const V SIUnits<V>::PROTON_MASS = 1.672621898e-27;

    template <typename V>
    const V SIUnits<V>::PROTON_CHARGE = SIUnits<V>::E;
  
    template <typename V>
    const V SIUnits<V>::roomT = 298.15;
  
    template <typename V>
    const V SIUnits<V>::kT = SIUnits<V>::KB * SIUnits<V>::roomT;
  
    template <typename V>
    const V SIUnits<V>::ALPHA = 7.2973525664e-03;
  
    template <typename V>
    const V SIUnits<V>::H = 6.626070040e-34;
  
    template <typename V>
    const V SIUnits<V>::MU_0 = 4.0 * MathConstants<V>::PI * 1.0e-07;
  
    template <typename V>
    const V SIUnits<V>::C_0 = 299792458;
  
    template <typename V>
    const V SIUnits<V>::eV = 
        std::sqrt(2.0 * SIUnits<V>::H * SIUnits<V>::ALPHA / (SIUnits<V>::MU_0 * SIUnits<V>::C_0));
  
    template <typename V>
    const V SIUnits<V>::WATER_VISCOSITY = 0.891*1.0e-03;
 
}

#endif /* SI_UNITS_HPP */

