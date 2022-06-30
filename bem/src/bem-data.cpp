/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 22:59.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/bem-data.hpp"

namespace simploce {

    BEMData::BEMData(const param_ptr_t& param, std::size_t nCol) :
        S{nCol, nCol}, lu{S}, b{nCol}, x{nCol}, epsRatio{}, ka{} {
        epsSolute = param->get<real_t>("bem.solute.eps");
        epsSolvent = param->get<real_t>("bem.solvent.eps");
        epsRatio = epsSolvent / epsSolute;
        ka = param->get<real_t>("bem.solvent.ka");
    }
    
}

