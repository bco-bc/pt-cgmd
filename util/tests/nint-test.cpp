/*
 * File:   nint-test.cpp
 * Author: André H. Juffer, Biocenter Oulu
 *
 * Created on 19 October 2019, 15:50
 */

#include "simploce/util/util.hpp"
#include "simploce/types/u-types.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

int main(int argc, char** argv) {
    std::cout << "nint-test test 1" << std::endl;
    real_t x0 = -1.55;
    real_t dx = 0.05;
    std::size_t N = 44;
    for (std::size_t k = 0; k != N; ++k) {
        real_t x = x0 + k * dx;
        std::cout << x << ' ' << util::nint(x) << std::endl;
    }
    std::cout << std::endl;

    real_t a = -9.5000000e-01;
    real_t n = util::nint(a);
    std::cout << a << ' ' << n << std::endl;
    return EXIT_SUCCESS;
}

