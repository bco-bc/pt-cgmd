/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 12:21.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_VOLTAGE_HPP
#define SIMULATION_VOLTAGE_HPP

#include "external-potential.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * An electric potential difference v(z) (voltage) over a distance z relative to a given reference point
     * at z=0 resulting in a constant homogeneous electric field E(z) = E0 in the z-direction. The
     * reference point is where the electric potential is zero, such that v(z) may be interpreted as the
     * electric potential at z.
     * For a positive voltage (that is, the electric potential is increasing in the positive z-direction relative to
     * the reference point), a constant electric field is pointing towards the negative z-axis. The force
     * on a positively (negatively) charged particle is in the same (opposite) direction. The electric
     * field E(r) = -grad(v(r)) = E0.  Since E0 is constant, one must have v(z)=-E0.z.
     */
    class Voltage : public external_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param voltage Voltage in the z-direction.
         * @param distance Distance for "voltage drop".
         * @param bc Boundary condition.
         * @param eps_r Relative permittivity, screens external field.
         * @param mesoscopic If true, this potential is for mesoscopic simulations (e.g., DPD).
         */
        Voltage(voltage_t voltage,
                dist_t distance,
                bc_ptr_t bc,
                real_t eps_r = 1.0,
                bool mesoscopic = false);

        /**
         * Constructor
         * @param E0 Component of the electric field.
         * @param eps_r Relative permittivity for screening.
         * @param bc Boundary condition.
         * @param mesoscopic mesoscopic If true, this potential is for mesoscopic simulations (e.g., DPD).
         */
        Voltage(el_field_t e0,
                bc_ptr_t bc,
                real_t eps_r = 1.0,
                bool mesoscopic = false);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        el_field_t e0_;
        bc_ptr_t bc_;
        real_t eps_r_;
        bool mesoscopic_;
    };
}

#endif //SIMULATION_VOLTAGE_HPP
