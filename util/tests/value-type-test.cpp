/*
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on August 9, 2019, 1:37 PM
 */

#include "simploce/types/value_t.hpp"
#include <cstdlib>
#include <iostream>

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "value_type_test test 1" << std::endl;
    
    using result_t = simploce::value_t<double, 0>;
    
    using pressure_t = simploce::value_t<double, 1>;
    using volume_t = simploce::value_t<double, 2>;
    using temperature_t = simploce::value_t<double, 3>;
    
    pressure_t p1{2};
    pressure_t p3 = p1;
    pressure_t p2(4);
    pressure_t p4{p2};
    std::cout << "p1: " << p1 << std::endl;
    std::cout << "p2: " << p2 << std::endl;
    std::cout << "p3 = p1: " << p3 << std::endl;
    std::cout << "p4{p2}: " << p4 << std::endl;
    std::cout << "(p2 == p4)? " << (p2==p4) << std::endl;
    std::cout << "(p2 == p3)? " << (p2==p3) << std::endl;
    pressure_t p = p1 + p2;
    std::cout << "p : " << p << std::endl;
        
    volume_t V(300);
    std::cout << "V : " << V << std::endl;
    pressure_t r1 = p * 2.0;
    std::cout << "p * 2.0: " << r1 << std::endl;    
    std::cout << "3.0 * V: " << 3.0 * V << std::endl;
    
    result_t r2 = p * V;
    std::cout << "p * V: " << r2 << std::endl;
    std::cout << "p()*V(): " << p() * V() << std::endl;
    
    std::cout << "2.0 * V * 3.0  : " << 2.0 * V * 3.0 << std::endl;
    std::cout << "2.0 * V() * 3.0: " << 2.0 * V() * 3.0 << std::endl;
    
    temperature_t T(300);
    
    std::cout << "p * V / T: " << p * V / T << std::endl;
    std::cout << "p() * V() / T(): " << p() * V() / T() << std::endl;
}

void test2() {
    using pressure_t = simploce::value_t<double, 1>;
    using temperature_t = simploce::value_t<double, 2>;

    pressure_t p{2.0};
    temperature_t T{4.0};
    double v = p * T;
    std::cout << "v = p * T = " << v << std::endl;
    std::cout << "v * (p * T) = " << v * (p * T) << std::endl;
}

void error() {
    using pressure_t = simploce::value_t<double, 1>;
    using volume_t = simploce::value_t<double, 2>;
    using result_t = simploce::value_t<double, 0>;

    pressure_t p{2.0};
    volume_t V{3.0};

    // Following will fail.
    //result_t r{p};
    //std::cout << "p + V = " << p + V << std::endl;
}


int main() {
    test1();
    test2();
    //error();
    return (EXIT_SUCCESS);
}

