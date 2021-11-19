/*
 * Author: André H. Juffer
 *
 * Created on 22 September 2019, 15:39
 */

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include "atypes.hpp"
#include "simploce/particle/bead.hpp"
#include <iostream>

namespace simploce {

    /**
     * Specifies an analysis of a given simulation model.
     */
    template <typename P>
    class Analysis;
    
    template <>
    class Analysis<Bead> {
    public:
        
        /**
         * Constructor
         * @param sm Simulation model
         * @param analyzer Analyzer.
         */
        Analysis(const cg_sim_model_ptr_t& sm,
                 const cg_analyzer_ptr_t& analyzer);
        
        /**
         * Drives the analysis.
         * @param param Parameters, such as the number of states in the trajectory
         * to skip (parameter name is "skip") before carrying out the analysis.
         * @param trajectory Input trajectory stream.
         */
        void perform(const sim_param_t& param,
                     std::istream& trajectory);
        
    private:
        
        cg_sim_model_ptr_t sm_;
        cg_analyzer_ptr_t analyzer_;
    };    
    
}

#endif /* ANALYSIS_HPP */

