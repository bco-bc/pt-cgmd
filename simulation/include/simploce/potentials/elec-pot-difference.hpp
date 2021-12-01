/*
 * Author: André H. Juffer.
 * Created on 30/11/2021, 12:21.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_ELEC_POT_DIFFERENCE_HPP
#define SIMULATION_ELEC_POT_DIFFERENCE_HPP

#include "external-potential.hpp"
#include "const-surface-charge-density.hpp"

namespace simploce {

    /**
     * An electric potential difference between two "points".
     */
    class ElectricPotentialDifference : public external_potential {
    public:

        /**
         * Direction along which the potential difference is applied.
         * x : Electric potential difference is applied along the x-direction.
         * y : Electric potential difference is applied along the y-direction.
         * z : Electric potential difference is applied along the z-direction.
         */
        enum DIRECTION {x = 1, y, z};

        /**
         * Conversion from char to DIRECTION.
         * @param value One of 'x', 'y', and 'z'.
         * @return DIRECTION.
         */
        static DIRECTION valueOf(char value);

        /**
         * Constructor. All arguments are required.
         * @param deltaV Electric potential difference. If positive and in the x-direction, an electric field is
         * pointing towards the negative x-axis. The force on a positively (negatively) charged particle is in
         * the same (opposite) direction.
         * @param distance Distance between the two points.
         * @param eps_r Relative permittivity.
         * @param bc Boundary condition.
         * @param direction Direction.
         */
        ElectricPotentialDifference(el_pot_diff deltaV,
                                    dist_t distance,
                                    real_t eps_r,
                                    bc_ptr_t bc,
                                    DIRECTION direction);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        bc_ptr_t bc_;
    };
}

#endif //SIMULATION_ELEC_POT_DIFFERENCE_HPP
