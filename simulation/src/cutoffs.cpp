/*
 * Author: André H. Juffer.
 * Created on 24/08/2023, 13:46
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/potentials/cutoffs.hpp"
#include "simploce/util/logger.hpp"

namespace simploce {

    Cutoffs::Cutoffs(dist_t sr, dist_t lr) :
        sr_{sr}, lr_{lr} {
            util::Logger logger("simploce::Cutoffs::Cutoffs(dist_t sr, dist_t lr)");
            if (sr_() <= 0.0 || lr_() <= 0.0) {
                std::string msg = "A cutoff distance must be > 0.0.";
                logAndThrow(logger, msg);
            }
    }

    dist_t
    Cutoffs::shortRanged() const {
        return sr_;
    }

    dist_t
    Cutoffs::longRanged() const {
        return lr_;
    }
}
