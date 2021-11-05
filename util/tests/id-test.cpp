//
// Created by ajuffer on 11/2/21.
//

#include "simploce/types/u-types.hpp"
#include "simploce/util/util.hpp"
#include <cstdlib>

int main() {
    simploce::id_t id = simploce::util::generateId();
    std::cout << "ID: " << simploce::util::toString(id) << std::endl;
    id = simploce::util::generateId();
    std::cout << "ID: " << simploce::util::toString(id) << std::endl;
    id = simploce::util::generateId();
    std::cout << "ID: " << simploce::util::toString(id) << std::endl;
    return (EXIT_SUCCESS);
}
