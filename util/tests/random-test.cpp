/*
 * Author: André H. Juffer.
 * Created on 21/06/22, 14:03.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/util.hpp"
#include <cstdlib>

int main() {
    for (auto i = 0; i != 500; ++i) {
        std::cout << simploce::util::randomUniform(0.0, 40.0) << std::endl;
    }
    return EXIT_SUCCESS;
}
