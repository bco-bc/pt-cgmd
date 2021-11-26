#include "simploce/analysis/analysis.hpp"
#include "simploce/analysis/analyzer.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"
#include <iostream>
#include <utility>

namespace simploce {
    
    Analysis::Analysis(p_system_ptr_t particleSystem,
                       a_param_ptr_t analysisParameters,
                       a_ptr_t analyzer) :
        particleSystem_{std::move(particleSystem)}, analysisParameters_{std::move(analysisParameters)},
        analyzer_{std::move(analyzer)} {
    }
        
    void 
    Analysis::perform(std::istream& trajectory)
    {
        util::Logger logger("simploce::Analysis::perform()");
        const auto nskip = analysisParameters_->get<std::size_t>("analysis.trajectory.nskip");
        logger.info("Skipping first " + util::toString(nskip) + " states in trajectory.");
        
        std::size_t counter = 0;
        particleSystem_->readState(trajectory);
        while (trajectory.good() ) {
            counter += 1;
            if ( counter > nskip ) {
                analyzer_->perform(particleSystem_);
            }
            particleSystem_->readState(trajectory);
        }
    }
    
}