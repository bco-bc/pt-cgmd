/*
 * File:   mu-units-test.cpp
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on October 31, 2019, 4:07 PM
 */

#include "simploce/units/units-mu.hpp"
#include "simploce/conf/u-conf.hpp"
#include <iostream>
#include <iomanip>

using namespace simploce;

int main(int argc, char** argv) {
    using V = double;
    std::cout.setf(std::ios::scientific);
    std::cout.precision(conf::PRECISION);
    std::cout << "Factors expressed in molecular units (mu):" << std::endl;
    std::cout << "E =                        " << std::setw(conf::REAL_WIDTH) << units::mu<V>::E << std::endl;
    std::cout << "E0 =                       " << std::setw(conf::REAL_WIDTH) << units::mu<V>::E0 << std::endl;
    std::cout << "4 * PI * E0 =              " << std::setw(conf::REAL_WIDTH) << units::mu<V>::FOUR_PI_E0 << std::endl;
    std::cout << "F_EL =                     " << std::setw(conf::REAL_WIDTH) << units::mu<V>::F_EL << std::endl;
    std::cout << "KB =                       " << std::setw(conf::REAL_WIDTH) << units::mu<V>::KB << std::endl;
    std::cout << "R =                        " << std::setw(conf::REAL_WIDTH) << units::mu<V>::R << std::endl;
    std::cout << "F =                        " << std::setw(conf::REAL_WIDTH) << units::mu<V>::F << std::endl;
    std::cout << "kT =                       " << std::setw(conf::REAL_WIDTH) << units::mu<V>::kT << std::endl;
    std::cout << "PROTON_MASS =              " << std::setw(conf::REAL_WIDTH) << units::mu<V>::PROTON_MASS << std::endl;
    std::cout << "PROTON_CHARGE =            " << std::setw(conf::REAL_WIDTH) << units::mu<V>::PROTON_CHARGE << std::endl;
    std::cout << "WATER_VISCOSITY =          " << std::setw(conf::REAL_WIDTH) << units::mu<V>::WATER_VISCOSITY << std::endl;
    std::cout << "l_to_nm3 =                 " << std::setw(conf::REAL_WIDTH) << units::mu<V>::l_to_nm3 << std::endl;
    std::cout << "Angstrom_to_nm =           " << std::setw(conf::REAL_WIDTH) << units::mu<V>::Angstrom_to_nm << std::endl;
    std::cout << "V_to_kJ_mol_e =            " << std::setw(conf::REAL_WIDTH) << units::mu<V>::V_to_kJ_mol_e << std::endl;
    std::cout << "e_nm_to_D =                " << std::setw(conf::REAL_WIDTH) << units::mu<V>::e_nm_to_D << std::endl;
    std::cout << "nm_to_m=                   " << std::setw(conf::REAL_WIDTH) << units::mu<V>::nm_to_m << std::endl;
    std::cout << "cal_to_J=                  " << std::setw(conf::REAL_WIDTH) << units::mu<V>::cal_to_J << std::endl;
    std::cout << "kcal_mol_A2_to_kJ_mol_nm2= " << std::setw(conf::REAL_WIDTH) << units::mu<V>::kcal_mol_A2_to_kJ_mol_nm2 << std::endl;
}

