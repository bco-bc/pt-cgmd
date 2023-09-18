/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 9/23/2022.
 */

#ifndef SIMULATION_PRESSURE_GRADIENT_HPP
#define SIMULATION_PRESSURE_GRADIENT_HPP

#include "external-potential.hpp"
#include "simploce/util/direction.hpp"

namespace simploce {

    /**
     * Pressure gradient is enforced by a constant force applied to each particle in a system.
     */
    class PressureGradient : public external_potential {
    public:

        /**
         * Constructor.
         * @param f Force applied.
         */
        explicit PressureGradient(force_t f);

        std::pair<energy_t, force_t> operator () (const p_ptr_t& particle) override;

    private:

        force_t f_;

    };
}

#endif //SIMULATION_PRESSURE_GRADIENT_HPP
