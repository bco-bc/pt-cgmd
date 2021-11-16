/*
 * Author: André H. Juffer.
 * Created on 10/11/2021, 19:56.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/param.hpp"
#include <cstdlib>

int main() {
    simploce::param::param_t param;

    param.put("debug.filename", "test.dat");
    param.add("debug.modules.module", "1.0");
    param.add("debug.modules.module", "2.0");

    simploce::param::write(std::cout, param);

    param.clear();

    param.put("simulation.nwrite", 10);
    param.put("simulation.npairlist", 5);
    simploce::param::write(std::cout, param);
}

