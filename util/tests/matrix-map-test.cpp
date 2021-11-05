/*
 * File:   matrix-map-test.cpp
 * Author: juffer
 *
 * Created on September 16, 2019, 8:05 PM
 */

#include "simploce/util/map2.hpp"
#include "simploce/util/map.hpp"
#include <iostream>

using namespace simploce;

using map2_t = MatrixMap<std::string, double>;

int main(int argc, char** argv) {
    std::cout << "matrix-map-test test 1" << std::endl;
    map2_t map2;
    map2.add("a1", "a2", 2.0);
    map2.add("a2", "a3", -567.90);
    std::cout << map2 << std::endl;

    std::cout << "a1, a2: " << map2.at("a1", "a2") << std::endl;

    // Error.
    try {
        std::cout << "Following should result in an exception..." << std::endl;
        std::cout << map2.at("a1", "a3") << std::endl;
    } catch(std::exception& ex) {
        std::cerr << ex.what() << std::endl;
    }
    auto keys = std::make_pair("test1", "test2");
    std::cout << map2.get(keys) << std::endl;

    std::map<std::string, double> mp;
    std::cout << mp << std::endl;

    return 0;
}

