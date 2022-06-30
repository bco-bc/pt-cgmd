/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 24, 2022, 12:09 PM
 */


#include "simploce/units/units-dpd.hpp"
#include <iostream>

using namespace simploce;

int main() {
    const auto ROOM_TEMPERATURE = units::si<real_t>::ROOM_TEMPERATURE;

    real_t nm3 = units::mu<real_t>::nm_to_m;
    nm3 *= (nm3 * nm3);

    // See Yiannourakou et al, 2012.
    mass_t mass{4.72e-25};
    length_t length{10.26e-10};
    energy_t energy{4.11e-21};
    real_t gamma{1.93e-13};  // kg/s
    real_t alpha{1.06e-10};

    std::cout << "Characteristic values in SI:" << std::endl;
    std::cout << "Mass: " << mass << " kg" << std::endl;
    std::cout << "Length: " << length << " m" << std::endl;
    std::cout << "Energy: " << energy << " J" << std::endl;
    std::cout << "Other values: " << std::endl;
    std::cout << "Friction coefficient (gamma): " << gamma << " kg/s" << std::endl;
    std::cout << "Maximum repulsion (i,j): : " << alpha << " J/m" << std::endl;
    std::cout << "kT: " << units::si<real_t>::kT << " J" << std::endl;
    std::cout << std::endl;

    mass /= units::si<real_t>::MU;
    length /= units::mu<real_t>::nm_to_m;
    energy = units::mu<real_t>::energyJtoKJperMol(energy());
    gamma /= units::si<real_t>::MU; // mu/s
    gamma *= 1.0e-12;
    alpha = units::mu<real_t>::energyJtoKJperMol(alpha);
    alpha *= 1.0e-09;
    position_t r{0.1, 0.2, 0.3};
    velocity_t v{ 1.0, 2.0, 3.0};
    number_density_t rho{2.78e+27 * nm3};
    stime_t dt{units::mu<real_t>::timeToPs(4.39e-13)};
    temperature_t temperature{ROOM_TEMPERATURE};

    std::cout << "Characteristic values in MU:" << std::endl;
    std::cout << "Mass: " << mass << " mu" << std::endl;
    std::cout << "Length: " << length << " nm" << std::endl;
    std::cout << "Energy: " << energy << " kJ/mol" << std::endl;
    std::cout << "Other values: " << std::endl;
    std::cout << "Time step: " << dt << " ps" << std::endl;
    std::cout << "Friction coefficient (gamma): " << gamma << " mu/ps" << std::endl;
    std::cout << "Maximum repulsion (i,j): : " << alpha << " kJ/(mol nm)" <<  std::endl;
    std::cout << "Position: " << r << " nm" << std::endl;
    std::cout << "Velocity: " << v << " nm/ps" << std::endl;
    std::cout << "Number density: " << rho << " nm⁻3" << std::endl;
    std::cout << "Room temperature: " << temperature << " K" << std::endl;
    std::cout << "kT: " << units::mu<real_t>::kT << " kJ/mol" << std::endl;
    std::cout << std::endl;

    units::dpd<real_t> dpd{mass, length,energy};

    std::cout << "MVV_DPD <-> MU units:" << std::endl;
    std::cout << "Mass (MVV_DPD): " << dpd.mass(mass) << std::endl;
    std::cout << "Mass (MU): " << dpd.inverseMass(dpd.mass(mass)) << " u" << std::endl;
    std::cout << "Length (MVV_DPD): " << dpd.length(length) << std::endl;
    std::cout << "Length (MU): " << dpd.inverseLength(dpd.length(length)) << " nm" << std::endl;
    std::cout << "Energy (MVV_DPD): " << dpd.energy(energy) << std::endl;
    std::cout << "Energy (MU): " << dpd.inverseEnergy(dpd.energy(energy)) << " kJ/mol" << std::endl;
    std::cout << "Number density (MVV_DPD): " << dpd.numberDensity(rho) << std::endl;
    std::cout << "Number density (MU): " << dpd.inverseNumberDensity(dpd.numberDensity(rho)) << std::endl;
    std::cout << "Time step (MVV_DPD): " << dpd.time(dt) << std::endl;
    std::cout << "Time step (MU): " << dpd.inverseTime(dpd.time(dt)) << " ps" << std::endl;
    std::cout << "Friction coefficient (MVV_DPD): " << dpd.gamma(gamma) << std::endl;
    std::cout << "Friction coefficient (MU): " << dpd.inverseGamma(dpd.gamma(gamma)) << " mu/s" << std::endl;
    std::cout << "Maximum repulsion (i,j) (MVV_DPD): " << dpd.alpha(alpha) << std::endl;
    std::cout << "Maximum repulsion (i,j) (MU): " << dpd.inverseAlpha(dpd.alpha(alpha)) << " kJ/(mol nm)" << std::endl;
    std::cout << "Position (MVV_DPD): " << dpd.position(r) << std::endl;
    std::cout << "Position (MU): " << dpd.inversePosition(dpd.position(r)) << " nm" << std::endl;
    std::cout << "Velocity (MVV_DPD): " << dpd.velocity(v) << std::endl;
    std::cout << "Velocity (MU): " << dpd.inverseVelocity(dpd.velocity(v)) << " nm/ps" << std::endl;
    std::cout << "Room temperature MVV_DPD: " << dpd.temperature(temperature) << std::endl;
    std::cout << "Room temperature (MU): " << dpd.inverseTemperature(dpd.temperature(temperature))
              << " K" << std::endl;

    temperature = 2.0 * dpd.temperature(temperature_t{ROOM_TEMPERATURE});
    std::cout << "Temperature MVV_DPD: " << temperature << std::endl;
    std::cout << "Temperature MU: " << dpd.inverseTemperature(temperature) << std::endl;

}