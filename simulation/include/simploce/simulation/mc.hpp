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
     * Monte Carlo (MC). Following parameters are applicable (JSON format, example values, with explanation):
     * {
     *  "simulation" : {
     *    "temperature" : "298.15",       Temperature
     *    "displacer": {
     *      "mc": {                       Displacer is Monte Carlo
     *          "range": "1.0",           Virtual box size around current position from which next position is selected.
     *          "keep-in-box": "x,y",     x- and y-coordinates are always in the box (if provided, ignore range for these coordinates).
     *          "in-box": "false",        Whether next particle position is randomly in the box (if true, ignore range for all x, y, z-coordinates).
     *          "z-non-negative": "true"  Whether the z-coordinate must be a non-negative number. If negative, perform reflection (z -> -z).
     *      }
     * }
     */
    class MonteCarlo : public displacer {
    public:
                
        /**
         * Constructor. All arguments are required.
         * @param params Simulation parameters.
         * @param interactor Interactor.
         */
        MonteCarlo(param_ptr_t param,
                   interactor_ptr_t interactor);

        /**
         * Displaces moveable particles, frozen particles are excluded.
         */
        SimulationData displace(const p_system_ptr_t& particles) const override;
        
    private:

        param_ptr_t param_;
        interactor_ptr_t interactor_;
    };
}

#endif /* MC_HPP */

