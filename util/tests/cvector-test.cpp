/*
 * File:   cvector-test.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 12:03 PM
 */

#include "simploce/types/cvector_t.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

int main(int argc, char** argv) {
    using fv_t = simploce::cvector_t<float, 2>;

    fv_t v1{1, 2, 3};
    fv_t v2{v1};
    std::cout << "v1: " << v1 << std::endl;
    std::cout << "v2: " << v1 << std::endl;
    std::cout << "|v1|: " << simploce::norm<double>(v1) << std::endl;
    std::cout << "|v1|²: " << simploce::norm_square<double>(v1) << std::endl;

    std::vector<float> v{2.10000, 2.20000, 2.30000};
    std:: cout << "v.size(): " << v.size() <<std::endl;
    fv_t vv{v};
    std::cout << "vv: " << vv << std::endl;

    std::array<double, 3> av{10.0, 11.0, -12.0};

    using position_t = simploce::cvector_t<double, 1>;
    position_t r1{1.0, 0, 0};
    position_t r2{0.0, 0.0, 2.0};
    position_t r3{av};

    std::cout << "r1: " << r1 << std::endl;
    std::cout << "r2: " << r2 << std::endl;
    std::cout << "r3: " << r3 << std::endl;
    r2.reset();
    std::cout << "r2 after reset: " << r2 << std::endl;
    std::cout << "av: ";
    for (int k = 0; k != 3; ++k) {
        std::cout << "  " << av[k];
    }
    std::cout << std::endl;
    std::cout << "r3: " << r3 << std::endl;
    av = r3.toArray();
    std::cout << "av (r3.toArray()): ";
    for (int k = 0; k != 3; ++k) {
        std::cout << "  " << av[k];
    }

    position_t rd = r2 - r1 + r3;
    std::cout << "rd: " << rd << std::endl;
    std::cout << "|rd| = " << simploce::norm<double>(rd) << std::endl;

    rd *= 4.0;
    std::cout << "rd: " << rd << std::endl;

    r1 = {1.0, 2.0, 3.0};
    r2 = {-1.0, -2.0, -3.0};
    auto a = simploce::angle<double>(r1, r2);
    std::cout << "Angle (r1, r2): " << a << std::endl;
    std::cout << "0.5 * PI      : " << 0.5 * std::acos(-1.L) << std::endl;

    auto cp = simploce::cross<double,2>(r1, r2);
    std::cout << "cp: " << cp << std::endl;
    if (simploce::norm<float>(cp) > 0.0 ) {
        std::cout << "Angle (r1, cp): " << simploce::angle<double>(r1, cp) << std::endl;
    }

    position_t p{3.9157429e+00, -7.7622789e-01, -4.5052286e+00};
    std::cout << "p: " << p << std::endl;
    auto p2 = simploce::inner<double>(p, p);
    std::cout << "inner(p, p): " << p2 << std::endl;
    p2 = 0;
    for (std::size_t k = 0; k != 3; ++k) {
        p2 += p[k] * p[k];
    }
    std::cout << "p*p: " << p2 << std::endl;

    return 0;
}

