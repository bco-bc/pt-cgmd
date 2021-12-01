/*
 * Author: André H. Juffer.
 * Created on 26/11/2021, 17:09.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/analysis/distance.hpp"
#include "simploce/simulation/bc.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-system.hpp"
#include <stdexcept>
#include <utility>

namespace simploce {

    d_ptr_t
    Distance::create(const a_param_ptr_t& analysisParameters, const bc_ptr_t& bc) {
        return std::make_shared<Distance>(analysisParameters, bc);
    }

    Distance::Distance(const a_param_ptr_t& analysisParameters, bc_ptr_t bc)
        : id_i_{}, id_j_{}, bc_{std::move(bc)} {
        try {
            id_i_ = analysisParameters->get<std::string>("analysis.distance.id_i");
            id_j_ = analysisParameters->get<std::string>("analysis.distance.id_j");
        } catch(std::exception exception) {
            throw std::domain_error("simploce::Distance(): Missing particle identifier(s).");
        }
    }

    void
    Distance::perform(const p_system_ptr_t& particleSystem) {
        static int counter = 0;

        counter += 1;
        auto pi = particleSystem->find(id_i_);
        auto pj = particleSystem->find(id_j_);
        if (pi == nullptr || pj == nullptr) {
            throw std::domain_error("(" + id_i_ + "," + id_j_ + "): no such particle(s).");
        }
        auto rij = bc_->apply(pi->position(), pj->position());
        auto Rij = norm<real_t>(rij);
        auto pair = std::make_pair(counter, dist_t(Rij));
        distances_.emplace_back(pair);
    }

    std::vector<std::pair<int, dist_t>>
    Distance::results() {
        return distances_;
    }

}

