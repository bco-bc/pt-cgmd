/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 22:59.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/bem-data.hpp"

namespace simploce {

    BEMData::BEMData(const param_ptr_t& param, std::size_t dimension) :
            S{dimension, dimension}, lu{S}, b{dimension}, x{dimension}, epsRatio{} {
        epsSolute = param->get<real_t>("bem.solute.eps");
        epsSolvent = param->get<real_t>("bem.solvent.eps");
        epsRatio = epsSolvent / epsSolute;
    }
    
}

