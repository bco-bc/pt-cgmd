/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on May 24, 2022, 12:09 PM
 */


#include "simploce/units/units-dpd.hpp"
#include "simploce/util/box.hpp"
#include "simploce/types/cvector_t.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include <iostream>

using namespace simploce;

void Yiannourakou() {
    const auto ROOM_TEMPERATURE = units::si<real_t>::ROOM_TEMPERATURE;

    real_t nm3 = units::mu<real_t>::nm_to_m;
    nm3 *= (nm3 * nm3);

    // See Yiannourakou et al, 2012.
    mass_t mass{4.72e-25};        // kg
    length_t length{10.26e-10};   // m
    energy_t energy{4.11e-21};    // J.

    real_t gamma{1.93e-13};  // kg/s
    real_t alpha{1.06e-10};  // J/m

    std::cout << "Yiannourakou et al, 2012" << std::endl;
    std::cout << "Characteristic values in SI:" << std::endl;
    std::cout << "Mass: " << mass << " kg" << std::endl;
    std::cout << "Length: " << length << " m" << std::endl;
    std::cout << "Energy: " << energy << " J" << std::endl;
    std::cout << "Other values: " << std::endl;
    std::cout << "Friction coefficient (gamma): " << gamma << " kg/s" << std::endl;
    std::cout << "Maximum repulsion (i,j): : " << alpha << " J/m" << std::endl;
    std::cout << "kT: " << units::si<real_t>::kT << " J" << std::endl;
    std::cout << std::endl;

    mass = units::mu<real_t>::mass_kg_to_u(mass());                // mu.
    length = units::mu<real_t>::length_m_to_nm(length());         // nm.
    energy = units::mu<real_t>::energy_J_to_kJ_per_mol(energy());  // kJ/mol
    gamma /= units::si<real_t>::MU; // mu/s
    gamma *= 1.0e-12;
    alpha = units::mu<real_t>::energy_J_to_kJ_per_mol(alpha);
    alpha *= 1.0e-09;
    position_t r{0.1, 0.2, 0.3};
    velocity_t v{ 1.0, 2.0, 3.0};
    number_density_t rho{2.78e+27 * nm3};
    stime_t dt{units::mu<real_t>::time_s_to_ps(4.39e-13)};
    temperature_t temperature{ROOM_TEMPERATURE};

    std::cout << "Characteristic values in MU:" << std::endl;
    std::cout << "Mass: " << mass << " mu" << std::endl;
    std::cout << "Length: " << length << " nm" << std::endl;
    std::cout << "Energy: " << energy << " kJ/mol" << std::endl;
    std::cout << "Other values: " << std::endl;
    std::cout << "Time step: " << dt << " ps" << std::endl;
    std::cout << "Time step: " << units::mu<real_t>::time_ps_to_s(dt()) << " s " << std::endl;
    std::cout << "Friction coefficient (gamma): " << gamma << " mu/ps" << std::endl;
    std::cout << "Maximum repulsion (i,j): : " << alpha << " kJ/(mol nm)" <<  std::endl;
    std::cout << "Position: " << r << " nm" << std::endl;
    std::cout << "Velocity: " << v << " nm/ps" << std::endl;
    std::cout << "Number density: " << rho << " nm⁻3" << std::endl;
    std::cout << "Room temperature: " << temperature << " K" << std::endl;
    std::cout << "kT: " << units::mu<real_t>::kT << " kJ/mol" << std::endl;
    std::cout << std::endl;

    units::dpd<real_t> dpd{mass, length,energy};

    std::cout << "DPD <-> MU, SI units:" << std::endl;
    std::cout << "Mass (DPD): " << dpd.mass(mass) << std::endl;
    mass = dpd.inverseMass(dpd.mass(mass));
    std::cout << "Mass (MU): " << mass << " u" << std::endl;
    std::cout << "Mass (SI): " << units::mu<real_t>::mass_u_to_kg(mass()) << " kg" << std::endl;
    std::cout << "Length (DPD): " << dpd.length(length) << std::endl;
    length = dpd.inverseLength(dpd.length(length));
    std::cout << "Length (MU): " << length << " nm" << std::endl;
    std::cout << "Length (SI): " << units::mu<real_t>::length_nm_to_m(length()) << " m" << std::endl;
    std::cout << "Energy (DPD): " << dpd.energy(energy) << std::endl;
    energy = dpd.inverseEnergy(dpd.energy(energy));
    std::cout << "Energy (MU): " << energy << " kJ/mol" << std::endl;
    std::cout << "Energy (SI): " << units::mu<real_t>::energy_kJ_per_mol_to_J(energy()) << " J" << std::endl;
    std::cout << "Number density (DPD): " << dpd.numberDensity(rho) << std::endl;
    rho = dpd.inverseNumberDensity(dpd.numberDensity(rho));
    std::cout << "Number density (MU): " << rho << std::endl;
    std::cout << "Number density (SI): " << units::mu<real_t>::number_density_nm3_m3(rho()) << " m^-3" << std::endl;
    std::cout << "Time step (DPD): " << dpd.time(dt) << std::endl;
    dt = dpd.inverseTime(dpd.time(dt));
    std::cout << "Time step (MU): " << dt << " ps" << std::endl;
    std::cout << "Time step (SI): " << units::mu<real_t>::time_ps_to_s(dt()) << " s" << std::endl;
    std::cout << "Friction coefficient (DPD): " << dpd.gamma(gamma) << std::endl;
    gamma = dpd.inverseGamma(dpd.gamma(gamma));
    std::cout << "Friction coefficient (MU): " << gamma << " mu/ps" << std::endl;
    gamma = units::mu<real_t>::mass_u_to_kg(gamma);
    gamma /= units::mu<real_t>::time_ps_to_s(1.0);
    std::cout << "Friction coefficient (SI): " << gamma << " kg/s" << std::endl;
    std::cout << "Maximum repulsion (i,j) (DPD): " << dpd.alpha(alpha) << std::endl;
    alpha = dpd.inverseAlpha(dpd.alpha(alpha));
    std::cout << "Maximum repulsion (i,j) (MU): " << alpha << " kJ/(mol nm)" << std::endl;
    alpha = units::mu<real_t>::energy_kJ_per_mol_to_J(alpha);
    alpha /= units::mu<real_t>::length_nm_to_m(1.0);
    std::cout << "Maximum repulsion (i,j) (MU): " << alpha << " J/m" << std::endl;
    std::cout << "Position (DPD): " << dpd.position(r) << std::endl;
    std::cout << "Position (MU): " << dpd.inversePosition(dpd.position(r)) << " nm" << std::endl;
    std::cout << "Velocity (DPD): " << dpd.velocity(v) << std::endl;
    std::cout << "Velocity (MU): " << dpd.inverseVelocity(dpd.velocity(v)) << " nm/ps" << std::endl;
    std::cout << "Room temperature (DPD): " << dpd.temperature(temperature) << std::endl;
    std::cout << "Room temperature (MU): " << dpd.inverseTemperature(dpd.temperature(temperature))
              << " K" << std::endl;

    temperature = 2.0 * dpd.temperature(temperature_t{ROOM_TEMPERATURE});
    std::cout << "Temperature (DPD): " << temperature << std::endl;
    std::cout << "Temperature (MU): " << dpd.inverseTemperature(temperature) << std::endl;
}

void polystyrene(const spec_catalog_ptr_t& catalog) {
    std::cout.setf(std::ios::scientific);

    real_t cm_to_m = 1.0e-02;
    real_t g_to_kg = 1.0e-03;
    auto H = catalog->H();
    auto C = catalog->C();

    std::cout << "Polystyrene PS:" << std::endl;
    std::cout << std::endl;
    std::cout << "SI units:" << std::endl;
    radius_t radiusBead = 3.0e-06;          // m.
    density_t rho = 1.04;                  // g/cm^3
    std::cout << "Radius single bead: " << radiusBead << " m" << std::endl;
    std::cout << "Density: " << rho << " g/cm^3" << std::endl;

    rho = rho() * g_to_kg;
    rho /= (cm_to_m * cm_to_m * cm_to_m);   // kg/m^3
    std::cout << "Density: " << rho << " kg/m^3" << std::endl;

    volume_t volumeBead = 4.0 * math::constants<real_t>::PI * radiusBead() * radiusBead() * radiusBead() / 3.0;
    std::cout << "Volume single bead: " << volumeBead << " m^3" << std::endl;
    mass_t massBead = volumeBead() * rho();
    std::cout << "Mass single bead: " << massBead << " kg" << std::endl;
    mass_t massBeadInU = massBead()/units::si<real_t>::MU;
    std::cout << "Mass single bead: " << massBeadInU << " mu" << std::endl;

    box_ptr_t box = factory::box(60e-06, 50e-06, 50e-03);
    std::cout << "Box size: " << *box << " m m m" << std::endl;
    std::cout << "Box volume: " << box->volume() << " m^3" << std::endl;

    mass_t massMonomerInU = 8.0 * H->mass() + 8.0 * C->mass();
    std::cout << "Mass one PS monomer: " << massMonomerInU << " u" << std::endl;
    real_t mMonomer = massMonomerInU();
    real_t mBead = massBeadInU();
    auto n = mBead / mMonomer;
    std::cout << "Number of monomers in single bead: " << size_t(n) << std::endl;
    real_t N_m = real_t(16.0 * n);
    std::cout << "Number of atoms in single bead: " << N_m << std::endl;

    std::cout << std::endl;
    std::cout << "DPD units:" << std::endl;
    real_t l = radiusBead(); // m.
    real_t m = massBead();  // mu
    real_t e = 0.0;
    std::cout << "Characteristic length l: " << l << " m" << std::endl;
    std::cout << "Characteristic mass m: " << m << " kg (" << massBead/units::si<real_t>::MU << " mu )" << std::endl;
    std::cout << "Characteristic energy e: " << e << " J" << std::endl;

    radiusBead /= l;
    std::cout << "Radius: " << radiusBead << std::endl;

    rho = rho() * (l * l * l) / m;
    std::cout << "Density: " << rho << std::endl;

    volumeBead /= (l * l * l);
    std::cout << "Volume single bead: " << volumeBead << std::endl;

    massBead /= m;
    std::cout << "Mass single bead: " << massBead << std::endl;

    std::cout << "Coarse-graining parameter: " << N_m << std::endl;

    box = factory::box(box->lengthX() / l, box->lengthY() / l, box->lengthZ() / l);
    std::cout << "Box size: " << *box << std::endl;
    std::cout << "Box volume: " << box->volume() << " m^3" << std::endl;
    std::cout << "Number of particles in box if number density=3.0: " << real_t(int(3.0 * box->volume())) << std::endl;

}

int main() {

    Yiannourakou();

    /*
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------------" << std::endl;
    auto catalog = factory::particleSpecCatalog("/localdisk/resources/particles-specs.dat");
    polystyrene(catalog);
    */
}