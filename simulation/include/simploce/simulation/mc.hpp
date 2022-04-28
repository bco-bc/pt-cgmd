/*
 * Author: Andre H. Juffer, Biocenter Oulu, University of Finland, Oulu.
 *
 * Created on 19 October 2019, 18:56
 */

#ifndef MC_HPP
#define MC_HPP

#include "displacer.hpp"
#include "simploce/types/s-types.hpp"
#include <iostream>

namespace simploce {
    
    /**
     * Monte Carlo (MC)
     */
    class MonteCarlo : public displacer {
    public:
                
        /**
         * Constructor. All arguments are required.
         * @param simulationParams Simulation parameters.
         * @param interactor Interactor.
         */
        MonteCarlo(param_ptr_t param,
                   interactor_ptr_t interactor);
        
        SimulationData displace(const p_system_ptr_t& particles) const override;
        
    private:

        param_ptr_t param_;
        interactor_ptr_t interactor_;
        
    };
}

#endif /* MC_HPP */

