/*
 * File:   units-mu.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 3:40 PM
 */

#ifndef MU_UNITS_HPP
#define MU_UNITS_HPP

#include "units-si.hpp"

namespace simploce {
    namespace units {

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
        template<typename V>
        struct mu {

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

            /**
             * e nm (dipole moment) to D (Debye)
             */
            static const V e_nm_to_D;

            /**
             * nm to m.
             */
            static const V nm_to_m;

            /**
             * Conversion factor from cal to J, according to thermodynamics.
             */
            static const V cal_to_J;

            /**
             * Conversion factor from kcal/(mol Å^2) to kJ/(mol nm^2).
             */
            static const V kcal_mol_A2_to_kJ_mol_nm2;

            /**
             * Converts energy in J to kJ/mol.
             * @param energy Energy in J.
             * @return energy in kJ/mol.
             */
            static V energyJtoKJperMol(V energy);

            /**
             * Converts time in s to time in ps.
             * @param time Time is s.
             * @return Time in ps.
             */
            static V timeToPs(V time);

        };

        template<typename V>
        const V mu<V>::E = 1.0;

        template<typename V>
        const V mu<V>::E0 =
                si<V>::E0 / (si<V>::E * 1.0e+09) * 1.0e+03 / (si<V>::eV * si<V>::NA);

        template<typename V>
        const V mu<V>::FOUR_PI_E0 = 4.0 * math::constants<V>::PI * mu<V>::E0;

        template<typename V>
        const V mu<V>::F_EL = 1.0 / mu<V>::FOUR_PI_E0;

        template<typename V>
        const V mu<V>::KB = si<V>::KB * si<V>::NA / 1.0e+03;

        template<typename V>
        const V mu<V>::R = si<V>::R / 1.0e+03;

        template<typename V>
        const V mu<V>::F = si<V>::F / si<V>::E;

        template<typename V>
        const V mu<V>::kT = mu<V>::KB * si<V>::ROOM_TEMPERATURE;

        template<typename V>
        const V mu<V>::PROTON_MASS = si<V>::PROTON_MASS / si<V>::MU;

        template<typename V>
        const V mu<V>::PROTON_CHARGE = 1.0;

        template<typename V>
        const V mu<V>::WATER_VISCOSITY =
                si<V>::WATER_VISCOSITY / (si<V>::MU * 1.0e+09 * 1.0e+12);

        template<typename V>
        const V mu<V>::l_to_nm3 = 1.0e-03 * 1.0e+27;

        template<typename V>
        const V mu<V>::Angstrom_to_nm = 0.1;

        template<typename V>
        const V mu<V>::V_to_kJ_mol_e = si<V>::NA * si<V>::E / 1.0e+03;

        template<typename V>
        const V mu<V>::nm_to_m = 1.0e-09;

        template<typename V>
        const V mu<V>::e_nm_to_D = si<V>::E * mu<V>::nm_to_m * 2.9979245817809e+29;

        template <typename V>
        const V mu<V>::cal_to_J = 4.184000;

        template <typename V>
        const V mu<V>::kcal_mol_A2_to_kJ_mol_nm2 = mu<V>::cal_to_J / (mu<V>::Angstrom_to_nm * mu<V>::Angstrom_to_nm);

        template<typename V>
        V mu<V>::energyJtoKJperMol(V energy) {
            return energy / 1000.0 * si<V>::NA;
        }

        template<typename V>
        V mu<V>::timeToPs(V time) {
            return time * 1.0e+12;
        }

    }
}

#endif /* MU_UNITS_HPP */

