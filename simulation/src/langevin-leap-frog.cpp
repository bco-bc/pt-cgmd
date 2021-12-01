/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on November 16, 2019, 5:39 PM
 */

#include "simploce/simulation/langevin-leap-frog.hpp"
#include <stdexcept>

namespace simploce {

    LangevinLeapFrog::LangevinLeapFrog(sim_param_ptr_t simulationParameters,
                                       interactor_ptr_t interactor) :
        simulationParameters_{std::move(simulationParameters)}, interactor_{std::move(interactor)}{
        if ( !simulationParameters_ ) {
            throw std::domain_error("LangevinLeapFrog: Missing simulation parameters.");
        }
        if ( !interactor_ ) {
            throw std::domain_error("LangevinLeapFrog: Missing interactor.");
        }
    }
        
    SimulationData 
    LangevinLeapFrog::displace(const p_system_ptr_t& particleSystem) const {
        throw std::domain_error("LangevinLeapFrogNot implemented yet.");
    }
    
        
}

