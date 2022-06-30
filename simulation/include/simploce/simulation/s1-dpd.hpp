/*
 * Author: André H. Juffer.
 * Created on 24/06/22, 13:32.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_S1_DPD_HPP
#define SIMULATION_S1_DPD_HPP

#include "displacer.hpp"

namespace simploce {

    /**
     * Dissipative particle dynamics. Implements the S1 (splitter) algorithm by Shardlow.
     * @see P. Nikunen et al., Computer Physics Communications 153 (2003) 407–423.
     */
    class S1_DPD : public displacer {
    public:

        /**
         * Constructor. Moves parameters.
         * @param param Simulation parameters.
         * @param interactor Interactor.
         * @param bc Boundary condition.
         * @param dpdUnits MVV_DPD Units convertor.
         */
        S1_DPD(param_ptr_t param,
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

#endif //SIMULATION_S1_DPD_HPP
