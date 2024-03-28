/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on 9/23/2022.
 */

#include "simploce/potentials/pressure-gradient.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/util.hpp"

namespace simploce {

    PressureGradient::PressureGradient(force_t f) : f_{f} {
    }

    std::pair<energy_t, force_t>
    PressureGradient::operator()(const simploce::p_ptr_t &particle) const {
        static util::Logger logger("simploce::PressureGradient::operator()");
        static int counter = 0;
        if (counter == 0)
            logger.info(util::to_string(f_) + ": External force due to pressure gradient.");
        auto r = particle->position();
        counter += 1;
        return std::make_pair(energy_t{-1.0 * inner<real_t>(f_, r)}, f_);
    }

}
