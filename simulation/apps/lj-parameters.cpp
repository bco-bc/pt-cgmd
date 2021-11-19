/*
 * Author: André H. Juffer.
 * Created on 18/11/2021, 19:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include <cstdlib>
#include "simploce/simulation/s-types.hpp"
#include "simploce/units/units-mu.hpp"
#include <iostream>

using namespace simploce;


void conversion1() {
    std::cout << "Rahman's Lennard-Jones interaction parameters 'sigma' and 'eps':" << std::endl;
    // From https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405"
    real_t sigma = 3.4;                            // Ångstrom.
    real_t eps = 120;                              // Kelvin.
    std::cout << "sigma = " << sigma << " Å" << std::endl;
    std::cout << "eps = " << eps << " K" << std::endl;

    sigma *= units::mu<real_t>::Angstrom_to_nm;    // nm.
    real_t sigma_6 = sigma * sigma * sigma * sigma * sigma * sigma;
    real_t sigma_12 = sigma_6 * sigma_6;
    eps *= units::mu<real_t>::KB;                 // kJ/mol.

    real_t C12 = 4.0 * eps * sigma_12;           // kJ nm^12/mol
    real_t C6 = 4.0 * eps * sigma_6;             // kJ nm^6/ mol

    std::cout << "Lennard-Jones interaction parameters C12 and C6:" << std::endl;
    std::cout << "C12 = " << C12 << " kJ nm^12/mol" << std::endl;
    std::cout << "C6 = " << C6 << " kJ nm^6/ mol" << std::endl;
}

int main() {

    conversion1();

    return EXIT_SUCCESS;
}
