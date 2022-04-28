/*
 * Author: André H. Juffer.
 * Created on 18/11/2021, 19:14.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include <cstdlib>
#include "simploce/types/s-types.hpp"
#include "simploce/units/units-mu.hpp"
#include <iostream>

using namespace simploce;

void sigmaEpsilonToC12C6(real_t sigma, real_t eps) {
    std::cout.precision(conf::PRECISION);
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
    std::cout << std::endl;
}


int main() {
    std::cout << "Rahman's Lennard-Jones interaction parameters 'sigma' and 'eps':" << std::endl;
    // From https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405"
    real_t sigma = 3.4;                            // Ångstrom.
    real_t eps = 120;                              // Kelvin.
    sigmaEpsilonToC12C6(sigma, eps);

    std::cout << "Lenard et al NaCl Lennard-Jones interaction parameters 'sigma' and 'eps':" << std::endl;
    // From Lenard et al, J. Chem. Phys., 126, 044509, 2007.
    std::cout << "Na+ - Na+" << std::endl;
    sigma = 2.443;                               // Ångstrom
    eps = 14.3288;                               // Kelvin
    sigmaEpsilonToC12C6(sigma, eps);
    std::cout << "Na+ - Cl-" << std::endl;
    sigma = 2.796;                               // Ångstrom
    eps = 42.4104;                               // Kelvin
    sigmaEpsilonToC12C6(sigma, eps);
    std::cout << "Cl- - Cl-" << std::endl;
    sigma = 3.487;                               // Ångstrom
    eps = 117.7604;                               // Kelvin
    sigmaEpsilonToC12C6(sigma, eps);


    return EXIT_SUCCESS;
}
