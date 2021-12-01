/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on November 16, 2019, 5:39 PM
 */

#ifndef LANGEVIN_LEAP_FROG_HPP
#define LANGEVIN_LEAP_FROG_HPP

#include "displacer.hpp"
#include "s-types.hpp"

namespace simploce {
    
    /**
     * Displaces particles according to a stochastic impulsive Langevin 
     * Leap-Frog for systems without constrains.
     * @see <a href="https://dx.doi.org/10.1021/ct3000876 | J.">
     *  Goga et al, J. Chem. Theory Comput. 2012, 8, 3637−3649.
     * </a>
     */
    class LangevinLeapFrog: public displacer {
    public:

        /**
         * Constructor. All arguments are required.
         * @param simulationParameters Simulation parameters.
         * @param interactor Interactor.
         */
        LangevinLeapFrog(sim_param_ptr_t simulationParameters,
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

#endif /* LANGEVIN_LEAP_FROG_HPP */

