/*
 * File:   pt-pairlist-test.cpp
 * Author: ajuffer
 *
 * Created on September 25, 2019, 3:54 PM
 */

#include "simploce/simulation/pt-pair-list-generator.hpp"
#include "simploce/simulation/sfactory.hpp"
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "pt-pairlist-test test 1" << std::endl;
    
    box_ptr_t box = std::make_shared<box_t>(7.0);
    bc_ptr_t bc = factory::pbc(box);
    ProtonTransferPairListGenerator g(bc);
    
    pt_pair_list_gen_ptr_t generator = 
            factory::protonTransferPairListGenerator(bc);
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% pt-pairlist-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (pt-pairlist-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (pt-pairlist-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

