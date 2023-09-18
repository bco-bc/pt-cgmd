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
     * An electric potential difference over a given distance in a given direction resulting in a constant electric
     * field in the given direction.
     */
    class Voltage : public external_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param voltage Electric potential difference. If positive and in the x/y/z-direction, a constant
         * electric field is pointing towards the negative x/y/z-axis. The force on a positively (negatively) charged particle is in
         * the same (opposite) direction. The electric field E0(r) = -grad(ep(r)) where ep(r) is the electric
         * potential at r and r is a position. Given that E0 being constant, p(r)=-E0.r (inner product)
         * @param distance Distance for "voltage drop".
         * @param eps_r Relative permittivity.
         * @param bc Boundary condition.
         * @param direction Direction.
         * @param mesoscopic If true, this potential is for mesoscopic simulations (e.g., DPD).
         */
        Voltage(voltage_t voltage,
                dist_t distance,
                real_t eps_r,
                const bc_ptr_t& bc,
                const Direction& direction,
                bool mesoscopic = false);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    };
}

#endif //SIMULATION_VOLTAGE_HPP
