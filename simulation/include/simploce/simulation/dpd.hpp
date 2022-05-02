/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on April 28, 2022, 15:11 PM
 */

#ifndef SIMULATION_DPD_HPP
#define SIMULATION_DPD_HPP

#include "displacer.hpp"

namespace simploce {

    /**
     * Dissipative particle dynamics. Employs the algorithm by
     * Groot and Warren, J. Chem. Phys., v. 107, p. 4423, 1997.
     */
    class DPD : public displacer {
    public:

        /**
         * Constructor.
         * @param param Simulation parameters.
         * @param interactor Interactor.
         * @param bc Boundary condition.
         */
        DPD(param_ptr_t param,
            interactor_ptr_t interactor,
            bc_ptr_t bc);

        SimulationData displace(const p_system_ptr_t& particleSystem) const override;

    private:

        param_ptr_t param_;
        interactor_ptr_t interactor_;
        bc_ptr_t bc_;
    };
}

#endif //SIMULATION_DPD_HPP
