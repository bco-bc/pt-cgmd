/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 20, 2019, 12:39 PM
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "simploce/types/s-types.hpp"
#include <fstream>

namespace simploce {
    
    /**
     * Performs a simulation. Writes particle state and simulation data to separate
     * output streams. Following parameters are applicable (JSON format, example values, with explanation):
     * {
     *   "simulation" : {
     *     "nsteps" : 1000,                            Total number of steps.
     *     "nwrite" : 10,                              Number of steps between writing data to output files (trajectory, simulation data)
     *     "temperature: : 298.15",                    Temperature.
     *     "mesoscale" : false,                        Whether mesoscale units are employed.
     *     "remove-com-motion" : false,                Whether the center-of-mass motion must be removed.
     *     "nremove-com-motion" "100",                 Number of steps between removal of center-of-mass motion, if applicable.
     *     "scale-velocities" : "false",               Whether velocities must be scaled to fit the temperature.
     *     "nscale-velocities" : 5,                    Number of steps between scaling velocities, if applicable.
     *     "relative-temperature-difference" : 10.0,   Relative difference or tolerance (in %) between requested and actual temperature, for velocity scaling.
     *     "freeze-boundary" : false,                  Whether surface boundary particles must be frozen (cannot move)
     *     "freeze: "SBP",                             Particle specification name of particles that must be frozen (cannot move).
     *     "bc" : {                                    Boundary condition (BC)
     *       "velocities" : {                          BC for velocities.
     *         "include" : false,                      Whether BC must be applied to particle velocities.
     *         "particle-groups" : false               Whether BC must be applied to velocities of particle at the particle group level.
     *       },
     *       "positions" : {                           BC for positions.
     *         "include" : false"                      Whether BC must be applied to particle positions.
     *       }
     *     }
     *   }
     * }
     */
    class Simulation {
    public:

        /**
         * Constructor. All arguments are required.
         * @param param Simulation parameters.
         * @param particleSystem Particular system.
         * @param catalog Particle specifications catalog.
         * @param displacer displacer.
         * @param bc Boundary condition.
         */
        explicit Simulation(param_ptr_t param,
                   p_system_ptr_t particleSystem,
                   spec_catalog_ptr_t catalog,
                   displacer_ptr_t displacer,
                   bc_ptr_t bc,
                   interactor_ptr_t interactor);
        
        /**
         * Performs the simulation.
         * @param trajectoryStream Output trajectory stream.
         * @param dataStream Output simulation data stream for monitoring the simulation.
         */
        void perform(std::ofstream& trajectoryStream,
                     std::ofstream& dataStream);
        
    private:
        
        param_ptr_t param_;
        p_system_ptr_t particleSystem_;
        spec_catalog_ptr_t catalog_;
        displacer_ptr_t displacer_;
        bc_ptr_t bc_;
        interactor_ptr_t interactor_;

    };


    
}

#endif /* SIMULATION_HPP */

