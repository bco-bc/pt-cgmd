/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on April 28, 2022, 15:11 PM
 */

#ifndef SIMULATION_MVV_DPD_HPP
#define SIMULATION_MVV_DPD_HPP

#include "displacer.hpp"

namespace simploce {

    /**
     * Dissipative particle dynamics. Implements the modified Velocity-Verlet (MVV) algorithm
     * by Groot and Warren, J. Chem. Phys., v. 107, p. 4423, 1997. It also employs a weight factor s
     * in (1-r/r_c)^s for calculating wR and wD. If not specified, s=1 (as in the original
     * DPD algorithm).
     * The units employed are MVV_DPD units, that is dimensionless (reduced) values are assumed.
     */
    class MVV_DPD : public displacer {
    public:

        /**
         * Constructor.
         * @param param Simulation parameters.
         * @param interactor Interactor.
         * @param bc Boundary condition.
         * @param dpdUnits DPD units convertor.
         */
        MVV_DPD(param_ptr_t param,
                interactor_ptr_t interactor,
                bc_ptr_t bc,
                units::dpd_ptr_t dpdUnits);

        SimulationData displace(const p_system_ptr_t& particleSystem) const override;

    private:

        param_ptr_t param_;
        interactor_ptr_t interactor_;
        bc_ptr_t bc_;
        units::dpd_ptr_t dpdUnits_;
    };
}

#endif //SIMULATION_MVV_DPD_HPP
