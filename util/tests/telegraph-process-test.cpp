/*
 * File:   telegraph-process-test.cpp
 * Author: ajuffer
 *
 * Created on October 4, 2019, 4:20 PM
 */

#include "simploce/util/telegraph-process.hpp"
#include <iostream>

using namespace simploce;

int main(int argc, char** argv) {
    auto results = TelegraphProcess::generate(1.0, 1.0, 100, 1);
    std::cout << results << std::endl;
    return 0;
}

