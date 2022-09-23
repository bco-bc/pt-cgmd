/*
 * Author: André H. Juffer.
 * Created on 21/06/22, 14:03.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/util.hpp"
#include <fstream>
#include <cstdlib>

using namespace simploce;

int main() {
    std::ofstream stream;
    std::string fileName = "/wrk3/tests/random.dat";
    stream.open(fileName, std::ios_base::out);
    for (auto i = 0; i != 2000000; ++i) {
        stream << simploce::util::randomUniform(0.0, 40.0) << std::endl;
    }
    stream.flush();
    stream.close();
    return EXIT_SUCCESS;
}
