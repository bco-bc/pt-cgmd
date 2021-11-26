/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 22 September 2019, 15:39
 */

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "a-types.hpp"
#include "../simulation/s-types.hpp"
#include "simploce/particle/bead.hpp"
#include <iostream>

namespace simploce {

    /**
     * Drives an analysis of a simulation. It reads states from a given trajectory and forwards them to the
     * analyzer.
     */
    class Analysis {
    public:
        
        /**
         * Constructor. All arguments are required.
         * @param particleSystem Particle system.
         * @param simulationParameters Analysis parameters. Must include the number of states in the trajectory
         * to skip (parameter name is "analysis.trajectory.nskip") before carrying out the analysis.
         * @param analyzer Analyzer.
         */
        Analysis(p_system_ptr_t particleSystem,
                 a_param_ptr_t analysisParameters,
                 a_ptr_t analyzer);
        
        /**
         * Drives the analysis.
         * @param trajectory Input trajectory stream.
         */
        void perform(std::istream& trajectory);
        
    private:

        a_param_ptr_t analysisParameters_;
        p_system_ptr_t particleSystem_;
        a_ptr_t analyzer_;
    };    
    
}

#endif /* ANALYSIS_HPP */

