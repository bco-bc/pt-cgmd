/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on August 16, 2019, 12:36 PM
 */

#ifndef LEAP_FROG_HPP
#define LEAP_FROG_HPP

#include "displacer.hpp"

namespace simploce {        
    
    /**
     * Standard leap frog algorithm for MD simulations.
     */
    class LeapFrog : public displacer {
    public:    
        
        LeapFrog(sim_param_ptr_t simulationParameters,
                 interactor_ptr_t interactor);
                
        /**
         * @return Time, potential energy (bonded, non-bonded), kinetic energy, temperature.
         */
        SimulationData displace(const p_system_ptr_t& particleSystem) const override;
        
    private:

        sim_param_ptr_t simulationParameters_;
        interactor_ptr_t interactor_;
    };
    

}

#endif /* LEAP_FROG_HPP */

