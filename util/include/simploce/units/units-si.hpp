/*
 * File:   units-si.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 3:34 PM
 */

#ifndef UNITS_SI_HPP
#define UNITS_SI_HPP

#include "simploce/util/math-constants.hpp"

namespace simploce {
    namespace units {

        /**
         * Numerical values of physical constants, in SI units (CODATA 2018).
         * @see <a href="http://physics.nist.gov/cuu/Constants/index.html">NIST</a>
         * @param V value type (float, double).
         */
        template<typename V>
        struct si {

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
             * Speed of light in vacuum. In m/s.
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

        // https://physics.nist.gov/cgi-bin/cuu/Value?e|search_for=elemtary+charge
        template<typename V>
        const V si<V>::E = 1.602176634e-19;

        // See https://physics.nist.gov/cgi-bin/cuu/Value?ep0|search_for=vacuum+permittivity
        template<typename V>
        const V si<V>::E0 = 8.8541878128e-12;

        template<typename V>
        const V si<V>::FOUR_PI_E0 = 4.0 * math::constants<V>::PI * si<V>::E0;

        template<typename V>
        const V si<V>::F_EL = 1.0 / si<V>::FOUR_PI_E0;

        // See https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=Boltzmann
        template<typename V>
        const V si<V>::KB = 1.380649e-23;

        // https://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=Avogadro
        template<typename V>
        const V si<V>::NA = 6.02214076e+23;

        template<typename V>
        const V si<V>::R = si<V>::KB * si<V>::NA;

        template<typename V>
        const V si<V>::F = si<V>::E * si<V>::NA;

        // See https://physics.nist.gov/cuu/Units/outside.html
        template<typename V>
        const V si<V>::MU = 1.66054e-27;

        // See https://physics.nist.gov/cgi-bin/cuu/Value?mp
        template<typename V>
        const V si<V>::PROTON_MASS = 1.67262192369e-27;

        template<typename V>
        const V si<V>::PROTON_CHARGE = si<V>::E;

        template<typename V>
        const V si<V>::roomT = 298.15;

        template<typename V>
        const V si<V>::kT = si<V>::KB * si<V>::roomT;

        // https://physics.nist.gov/cgi-bin/cuu/Value?alph|search_for=fine+structure+
        template<typename V>
        const V si<V>::ALPHA = 7.2973525693e-03;

        // See https://physics.nist.gov/cgi-bin/cuu/Value?h|search_for=Planck
        template<typename V>
        const V si<V>::H = 6.62607015e-34;

        template<typename V>
        const V si<V>::MU_0 = 4.0 * math::constants<V>::PI * 1.0e-07;

        // See https://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=Speed+of+light
        template<typename V>
        const V si<V>::C_0 = 299792458;

        template<typename V>
        const V si<V>::eV = std::sqrt(2.0 * si<V>::H * si<V>::ALPHA / (si<V>::MU_0 * si<V>::C_0));

        // See https://en.wikipedia.org/wiki/Properties_of_water
        template<typename V>
        const V si<V>::WATER_VISCOSITY = 0.890 * 1.0e-03;

    }
}

#endif /* UNITS_SI.HPP */

