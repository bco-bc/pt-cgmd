/*
 * Calculates total mass for several particle specifications.
 * File:  total-masses.cpp
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 18, 2019, 3:50 PM
 */

#include "simploce/util/file.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/p-types.hpp"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace simploce;

int main() {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);

    std::string fileName = "/localdisk/resources/particles-specs.dat";
    std::ifstream stream;
    util::open_input_file(stream, fileName);
    spec_catalog_ptr_t catalog = ParticleSpecCatalog::obtainFrom(stream);

    auto H = catalog->H();
    auto O = catalog->O();
    auto N = catalog->N();
    auto C = catalog->C();

    // H
    std::cout << "Mass H: " << H->mass() << std::endl;
    std::cout << "Mass O: " << O->mass() << std::endl;
    std::cout << "Mass N: " << N->mass() << std::endl;
    std::cout << "Mass C: " << C->mass() << std::endl;
    std::cout << std::endl;

    // pCOOH
    auto mass = C->mass() + 2.0 * O->mass();
    std::cout << "Mass specification pCOOH (deprotonated): " << mass << std::endl;
    std::cout << "Mass specification pCOOH (protonated):   " << mass + H->mass() << std::endl;
    std::cout << std::endl;

    // HCOOH
    mass = C->mass() + 2.0 * O->mass();
    std::cout << "Mass specification HCOOH (deprotonated): " << mass << std::endl;
    std::cout << "Mass specification HCOOH (protonated):   " << mass + H->mass() << std::endl;
    std::cout << std::endl;

    // mH2O
    mass = 2.0 * H->mass() + O->mass();
    std::cout << "Mass specification mH2O: " << mass << std::endl;
    std::cout << std::endl;

    // NH3
    mass = N->mass() + 3.0 * H->mass();
    std::cout << "Mass specification NH3: " << mass << std::endl;
    std::cout << std::endl;

    // NH4
    mass = N->mass() + 3.0 * H->mass();
    std::cout << "specification NH4 (deprotonated): " << mass << std::endl;
    std::cout << "specification NH4 (protonated):   " << mass + H->mass() << std::endl;
    std::cout << std::endl;

    return (EXIT_SUCCESS);
}