/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   telegraph-process-test.cpp
 * Author: ajuffer
 *
 * Created on October 4, 2019, 4:20 PM
 */

#include "simploce/util/telegraph-process.hpp"
#include <cstdlib>
#include <iostream>

using namespace simploce;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::clog << "telegraph-process-test test 1" << std::endl;
    
    auto results = TelegraphProcess::generate(1.0, 1.0, 100, 1);
    std::cout << results << std::endl;
}

void test2() {
    std::cout << "telegraph-process-test test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (telegraph-process-test) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% telegraph-process-test" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (telegraph-process-test)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (telegraph-process-test)" << std::endl;

    //std::cout << "%TEST_STARTED% test2 (telegraph-process-test)\n" << std::endl;
    //test2();
    //std::cout << "%TEST_FINISHED% time=0 test2 (telegraph-process-test)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

