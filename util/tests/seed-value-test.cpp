/*
 * File:   seed-value-test.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 21, 2019, 1:37 PM
 */

#include "simploce/util/util.hpp"
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
    for ( std::size_t counter = 0; counter != 10; ++counter) {
        auto value = simploce::util::seedValue<std::size_t>();
        std::cout << "Seed value: " << value << std::endl;
    }
    return (EXIT_SUCCESS);
}

