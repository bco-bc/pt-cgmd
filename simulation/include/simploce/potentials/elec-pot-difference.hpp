/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 12:21.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_ELEC_POT_DIFFERENCE_HPP
#define SIMULATION_ELEC_POT_DIFFERENCE_HPP

#include "external-potential.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * An electric potential difference over a given distance in a given direction
     */
    class ElectricPotentialDifference : public external_potential {
    public:

        /**
         * Constructor. All arguments are required.
         * @param deltaV Electric potential difference. If positive and in the x/y/z-direction, an electric field is
         * pointing towards the negative x/y/z-axis. The force on a positively (negatively) charged particle is in
         * the same (opposite) direction.
         * @param distance Distance.
         * @param eps_r Relative permittivity.
         * @param bc Boundary condition.
         * @param direction Direction.
         */
        ElectricPotentialDifference(el_pot_diff deltaV,
                                    dist_t distance,
                                    real_t eps_r,
                                    bc_ptr_t bc,
                                    const Direction& direction);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        real_t eps_r_;
        bc_ptr_t bc_;
    };
}

#endif //SIMULATION_ELEC_POT_DIFFERENCE_HPP
