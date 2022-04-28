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
        
        LeapFrog(param_ptr_t param,
                 interactor_ptr_t interactor);
                
        /**
         * @return Time, potential energy (bonded, non-bonded), kinetic energy, temperature.
         */
        SimulationData displace(const p_system_ptr_t& particles) const override;
        
    private:

        param_ptr_t param_;
        interactor_ptr_t interactor_;
    };
    

}

#endif /* LEAP_FROG_HPP */

