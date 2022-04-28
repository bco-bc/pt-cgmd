/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:38 PM
 */

#ifndef VELOCITY_VERLET_HPP
#define VELOCITY_VERLET_HPP

#include "displacer.hpp"
#include "sim-data.hpp"
#include "simploce/types/s-types.hpp"

namespace simploce {
    
    /**
     * Velocity Verlet algorithm for MD simulations.
     */
    class VelocityVerlet: public displacer {
    public:
        
        VelocityVerlet(param_ptr_t param,
                       interactor_ptr_t interactor);
        
        /**
         * @param at Atomistic model.
         * @return kinetic, potential energy, and temperature.
         */
        SimulationData displace(const p_system_ptr_t& particles) const override;
        
    private:
        
        param_ptr_t param_;
        interactor_ptr_t interactor_;
    };

}

#endif /* VELOCITY_VERLET_HPP */

