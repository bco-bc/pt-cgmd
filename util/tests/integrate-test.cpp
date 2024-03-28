//
// Created by juffer on 3/6/24.
//

#include "simploce/util/util.hpp"

using namespace simploce;

int main() {

    auto fn = [] (real_t x) {
        return x*x;
    };

    auto fn2 = [] (real_t x) {
        real_t t1 = 140000.0 / (140000.0 - 2100.0 * x);
        return 2000.0 * std::log(t1) -9.8 * x;
    };

    // First test, exact result is 2/3 = 0.66667.
    real_t result = util::integrate(fn, -1.0, 1.0);
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed);
    std::cout << std::setw(10) << result << " (exact: 2/3 = 0.66667)" << std::endl;

    // Second test. Exact result is 11061.34
    result = util::integrate(fn2, 8.0, 30.0);
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed);
    std::cout << std::setw(10) << result << " (exact: 11061.34)" << std::endl;
}
