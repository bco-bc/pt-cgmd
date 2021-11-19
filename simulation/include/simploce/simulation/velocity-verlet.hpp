/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 12:38 PM
 */

#ifndef VELOCITY_VERLET_HPP
#define VELOCITY_VERLET_HPP

#include "displacer.hpp"
#include "sim-data.hpp"
#include "s-types.hpp"

namespace simploce {
    
    /**
     * Velocity Verlet algorithm for MD simulations.
     */
    class VelocityVerlet: public Displacer {
    public:
        
        VelocityVerlet(sim_param_ptr_t simulationParameters,
                       interactor_ptr_t interactor);
        
        /**
         * Displaces beads of an atomistic model.
         * @param at Atomistic model.
         * @return kinetic, potential energy, and temperature.
         */
        SimulationData displace(const p_system_ptr_t& particleSystem) const override;
        
    private:
        
        sim_param_ptr_t simulationParameters_;
        interactor_ptr_t interactor_;
    };

}

#endif /* VELOCITY_VERLET_HPP */

