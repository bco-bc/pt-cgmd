/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 24, 2022, 11:41 AM
 */

#ifndef SIMULATION_UNITS_DPD_HPP
#define SIMULATION_UNITS_DPD_HPP

#include "simploce/types/u-types.hpp"
#include "simploce/units/units-mu.hpp"
#include <cmath>

namespace simploce {
    namespace units {

        /**
         * Converts molecular units (MU) to/from MVV_DPD (reduced, dimensionless) units, using characteristic values for
         * mass, length, and energy. All MVV_DPD quantities are expressed in these values, typically symbolically
         * referred to as 'm', 'l', and 'e' for units of mass, length, and energy, respectively. Time is expressed
         * in 't'. Temperature is expressed as kT=n e, where n is an non-negative integer and k is the Boltzmann
         * constant.
         * @see <a href="https://iopscience.iop.org/article/10.1209/0295-5075/97/34007/">Yiannourakou2012</a>
         * @tparam V Real value type.
         */
        template <typename V>
        class dpd {
        public:

            /**
             * Constructor
             * @param mass Characteristic mass in MU.
             * @param length Characteristic length in MU.
             * @param energy Characteristic energy in MU.
             */
            dpd(mass_t mass, length_t length, energy_t energy);

            /**
             * Returns mass in MVV_DPD units.
             * @param mass in MU.
             * @return Mass in MVV_DPD units.
             */
            mass_t mass(const mass_t& mass);

            /**
             * Returns mass in MU.
             * @param mass Mass in MVV_DPD units.
             * @return Mass in MU.
             */
            mass_t inverseMass(const mass_t& mass);

            /**
             * Returns length in MVV_DPD units.
             * @param length Length in MU.
             * @return Length in MVV_DPD units.
             */
            length_t length(const length_t& length);

            /**
             * Returns length in MU.
             * @param length Length in MVV_DPD units.
             * @return Length in MU
             */
            length_t inverseLength(const length_t& length);

            /**
             * Returns energy in MVV_DPD units.
             * @param energy Energy in MU.
             * @return Energy in MVV_DPD units.
             */
            energy_t energy(const energy_t& energy);

            /**
             * Returns energy in MU.
             * @param energy Energy in MVV_DPD units.
             * @return Energy in MU.
             */
            energy_t inverseEnergy(const energy_t& energy);

            /**
             * Returns time in MVV_DPD units.
             * @param time Time in MU.
             * @return Time in t
             */
            stime_t time(const stime_t& time);

            /**
             * Returns time in MU.
             * @param time Time in MVV_DPD units.
             * @return Time in MU.
             */
            stime_t inverseTime(const stime_t& time);

            /**
             * Returns maximum repulsion parameter in MVV_DPD units.
             * @param alpha Maximum repulsion parameter in MU.
             * @return Maximum repulsion parameter in MVV_DPD units.
             */
            V alpha(V alpha);

            /**
             * Returns maximum repulsion parameter in MU.
             * @param alpha Maximum repulsion parameter in MVV_DPD units.
             * @return Maximum repulsion parameter in MU.
             */
            V inverseAlpha(V alpha);

            /**
             * Returns friction coefficient in MVV_DPD units.
             * @param gamma Friction coefficient.
             * @return Friction coefficient in MVV_DPD units.
             */
            V gamma(V gamma);

            /**
             * Returns friction coefficient in MU.
             * @param gamma Friction coefficient in MVV_DPD.
             * @return Friction coefficient in MU.
             */
            V inverseGamma(V gamma);

            /**
             * Returns number density in MVV_DPD units.
             * @param rho Number density in MU units.
             * @return Number density in MVV_DPD units.
             */
            number_density_t numberDensity(const number_density_t& rho);

            /**
             * Returns number density in MU.
             * @param rho Number density in MVV_DPD units.
             * @return Number density in MU.
             */
            number_density_t inverseNumberDensity(const number_density_t& rho);

            /**
             * Returns temperature in MVV_DPD units.
             * @param temperature Temperature in MU.
             * @return Temperature in MVV_DPD units.
             */
            temperature_t temperature(const temperature_t& temperature);

            /**
             * Returns temperature in MU.
             * @param temperature Temperature in MVV_DPD units.
             * @return Temperature in MU.
             */
            temperature_t inverseTemperature(const temperature_t& temperature);

            /**
             * Returns position in MVV_DPD units.
             * @param r Position in MU.
             * @return Position in MVV_DPD units.
             */
            position_t position(const position_t& r);

            /**
             * Returns position in MU.
             * @param r Position in MVV_DPD units.
             * @return Position in MU.
             */
            position_t inversePosition(const position_t& r);

            /**
             * Returns velocity in MVV_DPD units.
             * @param v Velocity in MU.
             * @return Velocity in MVV_DPD units.
             */
            velocity_t velocity(const velocity_t& v);

            /**
             * Returns velocity in MU.
             * @param v Velocity in MVV_DPD units.
             * @return Velocity in MU.
             */
            velocity_t inverseVelocity(const velocity_t& v);

        private:

            mass_t mass_;
            length_t length_;
            energy_t energy_;

        };

        template <typename V>
        dpd<V>::dpd(mass_t mass, length_t length, energy_t energy) :
        mass_{mass}, length_(length), energy_{energy} {
        }

        template <typename V>
        mass_t
        dpd<V>::mass(const mass_t& mass) {
            static auto cf = 1.0 / mass_();
            return cf * mass();
        }

        template <typename V>
        mass_t
        dpd<V>::inverseMass(const mass_t& mass) {
            static auto cf = mass_();
            return cf * mass;
        }

        template <typename V>
        length_t
        dpd<V>::length(const length_t &length) {
            static auto cf = 1.0 / length_();
            return cf * length;
        }

        template <typename V>
        length_t
        dpd<V>::inverseLength(const length_t& length) {
            static auto cf = length_();
            return cf * length;
        }

        template <typename V>
        energy_t
        dpd<V>::energy(const energy_t &energy) {
            static auto cf = 1.0 / energy_();
            return cf * energy;
        }

        template <typename V>
        energy_t
        dpd<V>::inverseEnergy(const energy_t& energy) {
            static auto cf = energy_();
            return cf * energy;
        }

        template <typename V>
        stime_t
        dpd<V>::time(const stime_t &time) {
            // Follows from F=ma, where F has units J/m (kJ/mol) and ma has units kg m/s^2 (u nm/ps²).
            static auto cf = 1.0 / (length_() *  std::sqrt(mass_() / energy_()));
            return cf * time;
        }

        template <typename V>
        stime_t
        dpd<V>::inverseTime(const stime_t& time) {
            static auto cf = length_() *  std::sqrt(mass_() / energy_());
            return cf * time;
        }

        template<typename V>
        V
        dpd<V>::alpha(V alpha) {
            static auto cf = 1.0 / (energy_() / length_());
            return cf * alpha;
        }

        template <typename V>
        V
        dpd<V>::inverseAlpha(V alpha) {
            static auto cf = energy_() / length_();
            return cf * alpha;
        }

        template <typename V>
        V
        dpd<V>::gamma(V gamma) {
            static auto cf = 1.0 / (std::sqrt(mass_() * energy_()) / length_());
            return cf * gamma;
        }

        template <typename V>
        V
        dpd<V>::inverseGamma(V gamma) {
            static auto cf = std::sqrt(mass_() * energy_()) / length_();
            return cf * gamma;
        }

        template <typename V>
        number_density_t
        dpd<V>::numberDensity(const number_density_t &rho) {
            static auto cf = length_() * length_() * length_();
            return cf * rho;
        }

        template <typename V>
        number_density_t
        dpd<V>::inverseNumberDensity(const number_density_t& rho) {
            static auto cf = 1.0 / (length_() * length_() * length_());
            return cf * rho;
        }

        template <typename V>
        temperature_t
        dpd<V>::temperature(const temperature_t &temperature) {
            static auto cf = units::mu<V>::KB / energy_();
            return cf * temperature;
        }

        template <typename V>
        temperature_t
        dpd<V>::inverseTemperature(const temperature_t& temperature) {
            static auto cf = 1.0 / (units::mu<V>::KB / energy_());
            return cf * temperature;
        }

        template <typename V>
        position_t
        dpd<V>::position(const position_t& r) {
            static auto cf = 1.0 / length_();
            return cf * r;
        }

        template <typename V>
        position_t
        dpd<V>::inversePosition(const position_t& r) {
            static auto cf = length_();
            return cf * r;
        }

        template <typename V>
        velocity_t
        dpd<V>::velocity(const velocity_t& v) {
            static auto cf = this->time(1.0)() / length_();
            return cf * v;
        }

        template <typename V>
        velocity_t
        dpd<V>::inverseVelocity(const velocity_t& v) {
            static auto cf = 1.0 / (this->time(1.0)() / length_());
            return cf * v;
        }

    }
}

#endif //SIMULATION_UNITS_DPD_HPP
