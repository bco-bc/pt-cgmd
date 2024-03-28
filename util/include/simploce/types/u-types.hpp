/*
 * File:   u-types.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 4, 2019, 3:38 PM
 */

#ifndef U_TYPES_HPP
#define U_TYPES_HPP

#include "value_t.hpp"
#include "cvector_t.hpp"
#include "boost/uuid/uuid.hpp"
#include <memory>

namespace simploce {

    /**
     * Real numbers.
     */
    using real_t = double;

    /**
     * Identifier type.
     */
    //using id_t = boost::uuids::uuid;
    using id_t = std::string;

    /**
     * Energy.
     */
    using energy_t = value_t<real_t, -1>;

    /**
     * Time.
     */
    using stime_t = value_t<real_t, -2>;

    /**
     * Rate.
     */
    using rate_t = value_t<real_t, -3>;

    /**
     * Charge.
     */
    using charge_t = value_t<real_t, -4>;

    /**
     * pKa.
     */
    using pKa_t = value_t<real_t, -5>;

    /**
     * Mass.
     */
    using mass_t = value_t<real_t, -7>;

    /**
     * Density type (mass per unit volume).
     */
    using density_t = value_t<real_t, -8>;

    /**
     * Temperature.
     */
    using temperature_t = value_t<real_t, -9>;

    /**
     * Length or distance.
     */
    using length_t = value_t<real_t, -10>;
    using dist_t = length_t;

    /*
     * Radius and diameter.
     */
    using radius_t = length_t;
    using diameter_t = length_t;
    /**
     * Volume.
     */
    using volume_t = value_t<real_t, -11>;

    /**
     * Surface charge density (e.g. e/nm^2),
     */
    using srf_charge_density_t = value_t<real_t, -12>;

    /**
     * Relative permittivity.
     */
    using rel_perm_t = value_t<real_t, -13>;

    /**
     * Pressure.
     */
    using pressure_t = value_t<real_t, -14>;

    /**
     * Number density.
     */
    using number_density_t = value_t<real_t, -15>;

    /**
     * Molarity (amount of substance per unit volume, e.g. mol/L).
     */
    using molarity_t = value_t<real_t, -16>;

    /**
     * Molality is the number of moles of solute in a solution corresponding to
     * 1 kg or 1000 g of solvent.
     */
    using molality_t = value_t<real_t, -17>;

    /**
     * Electric potential.
     */
    using el_pot_t = value_t<real_t, -18>;

    /**
     * Electric potential difference.
     */
    using el_pot_diff_t = el_pot_t;
    using voltage_t = el_pot_diff_t;

    /**
     * Area
     */
     using area_t = value_t<real_t, -19>;

    /**
     * Position.
     */
    using position_t = cvector_t<real_t, 1>;
    using position_ptr_t = std::shared_ptr<position_t>;

    /*
     * Distance vector type.
     */
    using dist_vect_t = position_t;

    /**
     * Velocity.
     */
    using velocity_t = cvector_t<real_t, 2>;

    /**
     * Momentum.
     */
    using momentum_t = cvector_t<real_t, 3>;

    /**
     * Force.
     */
    using force_t = cvector_t<real_t, 4>;

    /**
     * Electric dipole moment.
     */
    using dipole_moment_t = cvector_t<real_t, 5>;

    /**
     * Electric field.
     */
    using el_field_t = cvector_t<real_t, 6>;

    /**
     * Normal vector.
     */
    using normal_t = dist_vect_t;

}

#endif /* U_TYPES_HPP */

