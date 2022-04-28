/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 16, 2019, 4:04 PM
 */

#ifndef LANGEVIN_VELOCITY_VERLET_HPP
#define LANGEVIN_VELOCITY_VERLET_HPP

#include "displacer.hpp"

namespace simploce {
    
    
    /**
     * Displaces particles according a stochastic velocity Verlet algorithm applicable to
     * an Langevin equation. Provides for a canonical ensemble (NVT constant) simulation.
     * @see <a href="http://dx.doi.org/10.1080/00268976.2012.760055">
     *   Grønbech-Jensen and Oded Farago, Molec. Phys., 111, 983-991, 2013
     * </a>
    */
    class LangevinVelocityVerlet : public displacer {
    public:

        /**
         * Constructor. All arguments are required.
         * @param param Simulation parameters.
         * @param interactor Interactor.
         */
        LangevinVelocityVerlet(param_ptr_t param,
                               interactor_ptr_t interactor);

        SimulationData displace(const p_system_ptr_t& particles) const override;

    private:

        param_ptr_t param_;
        interactor_ptr_t interactor_;
        
    };

}

#endif /* LANGEVIN_VELOCITY_VERLET_HPP */

