//
// Created by ajuffer on 5/18/22.
//

#include "simploce/util/util.hpp"
#include "simploce/types/u-types.hpp"
#include <iostream>

using namespace simploce;

int main() {

    real_t r1 = -2.0989;
    real_t r2 = 456.981;
    real_t r3 = 1.099;

    std::cout << "Sign r1: " << util::sgn(r1) << std::endl;
    std::cout << "Sign r2: " << util::sgn(r2) << std::endl;
    std::cout << "Sign r3: " << util::sgn(r3) << std::endl;

    auto same = util::sgn(r1)==util::sgn(r2);
    std::cout << "r1 and r2 have same sign? " << same << std::endl;
    same = util::sgn(r1)==util::sgn(r3);
    std::cout << "r1 and r3 have same sign? " << same << std::endl;
    same = util::sgn(r2)==util::sgn(r3);
    std::cout << "r2 and r3 have same sign? " << same << std::endl;
}
