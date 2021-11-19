/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 19 October 2019, 18:56
 */

#ifndef MC_HPP
#define MC_HPP

#include "displacer.hpp"
#include "s-types.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Monte Carlo (MC)
     */
    class MonteCarlo : public Displacer {
    public:
                
        /**
         * Constructor. All arguments are required.
         * @param simulationParams Simulation parameters.
         * @param interactor Interactor.
         */
        MonteCarlo(sim_param_ptr_t simulationParameters,
                   interactor_ptr_t interactor);
        
        SimulationData displace(const p_system_ptr_t& particleSystem) const override;
        
    private:

        sim_param_ptr_t simulationParameters_;
        interactor_ptr_t interactor_;
        
    };
}

#endif /* MC_HPP */

