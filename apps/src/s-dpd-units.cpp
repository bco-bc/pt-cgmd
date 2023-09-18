/*
 * Units for DPD simulations.
 */

#include "simploce/units/units-dpd.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/box.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/math-constants.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-spec.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <cmath>

namespace po = boost::program_options;
using namespace simploce;

namespace simploce {

    const real_t cm_to_m = 1.0e-02;
    const real_t g_to_kg = 1.0e-03;
    const real_t nm_to_m = units::mu<real_t>::nm_to_m;


    units::dpd<real_t>
    polystyreneInWaterChannel(const spec_catalog_ptr_t &catalog) {

        auto H = catalog->H();
        auto C = catalog->C();

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "POLYSTYRENE (PS) \"bead\":" << std::endl;
        std::cout << std::endl;
        std::cout << "SI units:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        length_t diameterBead = 6.0e-06;       // m.
        density_t rho = 1.04;                  // g/cm^3, from Wikipedia
        std::cout << "Diameter single bead (m):          " << diameterBead << std::endl;
        std::cout << "Density (g/cm^3):                  " << rho << std::endl;

        rho = rho() * g_to_kg;
        rho /= (cm_to_m * cm_to_m * cm_to_m);   // kg/m^3
        std::cout << "Density (kg/m^3):                  " << rho << std::endl;

        radius_t radiusBead = diameterBead() / 2.0;
        volume_t volumeBead = 4.0 * math::constants<real_t>::PI * radiusBead() * radiusBead() * radiusBead() / 3.0;
        std::cout << "Volume single bead (m^3):          " << volumeBead << std::endl;
        mass_t massBead = volumeBead() * rho();
        std::cout << "Mass single bead (kg):             " << massBead << std::endl;
        mass_t massBeadInU = units::mu<real_t>::mass_kg_to_u(massBead());
        std::cout << "Mass single bead (u):              " << massBeadInU << std::endl;

        mass_t massMonomerInU = 8.0 * H->mass() + 8.0 * C->mass();
        std::cout << "Mass single PS monomer (u):        " << massMonomerInU << std::endl;
        real_t mMonomer = massMonomerInU();
        real_t mBead = massBeadInU();
        auto n = mBead / mMonomer;
        std::cout << "Number of monomers in single bead:   " << std::round(n) << std::endl;
        std::cout << "Number of atoms per PS monomer:      " << 16 << std::endl;
        auto N_m = real_t(16.0 * n);
        std::cout << "Number of atoms in single bead:      " << std::round(N_m) << std::endl;
        std::cout << "Room temperature (K):                " << units::si<real_t>::ROOM_TEMPERATURE << std::endl;

        std::cout << std::endl;
        std::cout << "DPD units:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        real_t l = diameterBead(); // m.
        real_t m = massBead();     // mu
        real_t e = units::si<real_t>::kT;  // At room temperature, SATP.
        std::cout << "Characteristic length l:                          " << l << " m" << std::endl;
        std::cout << "Characteristic mass m:                            " << m << " kg (" << massBead / units::si<real_t>::MU << " mu)"
                  << std::endl;
        std::cout << "Characteristic energy e:                          " << e << " J (" << units::mu<real_t>::kT << " kJ/mol)" << std::endl;

        diameterBead /= l;
        std::cout << "Diameter:                                         " << diameterBead() << std::endl;

        rho = rho() * (l * l * l) / m;
        std::cout << "Density:                                          " << rho() << std::endl;

        volumeBead /= (l * l * l);
        std::cout << "Volume single bead:                               " << volumeBead() << std::endl;

        massBead /= m;
        std::cout << "Mass single bead:                                 " << massBead() << std::endl;

        std::cout << "Coarse-graining parameter:                        " << N_m << std::endl;

        std::cout << "Temperature kT=2 means T=" << 2.0 * e / units::si<real_t>::KB << " K" << std::endl;

        l = units::mu<real_t>::length_m_to_nm(l);
        m = units::mu<real_t>::mass_kg_to_u(m);
        e = units::mu<real_t>::energy_J_to_kJ_per_mol(e);
        return {m, l, e};
    }

    /**
     * Returns coarse-graining parameter.
     * @param dpd
     * @param catalog
     * @return Coarse-graining parameter.
     */
    real_t water(const units::dpd<real_t >& dpd, const spec_catalog_ptr_t &catalog) {
        auto water = catalog->molecularWater();
        mass_t massWaterMoleculeInU = water->mass();  // In u.
        mass_t massWaterMolecule = massWaterMoleculeInU * units::si<real_t>::MU;

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "WATER:" << std::endl;
        std::cout << "SI units:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "Using characteristic values:" << std::endl;
        std::cout << "Length (nm):                                  " << dpd.length()() << std::endl;
        std::cout << "Energy (kJ/mol):                              " << dpd.energy()() << std::endl;
        std::cout << "Mass (u):                                     " << dpd.mass()() << std::endl;

        density_t density{0.9970479};  // g/cm^3
        std::cout << "Water density at " << units::si<real_t>::ROOM_TEMPERATURE << " K (g/cm^3): " << density << std::endl;
        density = density() * (g_to_kg / (cm_to_m * cm_to_m * cm_to_m));
        std::cout << "Water density at " << units::si<real_t>::ROOM_TEMPERATURE << " K (kg/m^3): " << density << std::endl;
        std::cout << "Mass one water molecule (u):                  " << massWaterMoleculeInU() << std::endl;
        std::cout << "Mass one water molecule (kg):                 " << massWaterMolecule() << std::endl;

        length_t l = dpd.length() * nm_to_m;         // in m.
        radius_t radius = 0.5 * l();                 // in m.
        volume_t volume = 4.0 * math::constants<real_t>::PI * radius() * radius() * radius() / 3.0;  // in m^3
        std::cout << "Radius water bead (m):                        " << radius() << std::endl;
        std::cout << "Volume water bead (m^3):                      " << volume() << std::endl;

        mass_t massWaterBead = density() * volume();  // kg
        mass_t massWaterBeadInU = units::mu<real_t>::mass_kg_to_u(massWaterBead());
        std::cout << "Mass water bead (kg):                         " << massWaterBead() << std::endl;
        std::cout << "Mass water bead (u):                          " << massWaterBeadInU() << std::endl;
        number_density_t numberDensity = density/massWaterMolecule;
        std::cout << "Number density (m^-3):                        " << numberDensity() << std::endl;

        real_t numberOfWaterMolecules = massWaterBeadInU() / massWaterMoleculeInU();
        std::cout << "Number of water molecules in water bead:      " << std::size_t(numberOfWaterMolecules) << std::endl;

        std::cout << "DPD units:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;

        radius /= l();
        std::cout << "Diameter water bead:                          " << 2.0 * radius() << std::endl;

        massWaterBead = massWaterBeadInU / dpd.mass();
        std::cout << "Mass:                                         " << massWaterBead() << std::endl;
        std::cout << "Number of water molecules in water bead:      " << numberOfWaterMolecules << std::endl;
        number_density_t numberDensityDPD = 3.0;  // Desired number density.
        auto d = dpd.length() * 1.0e-09; // In m.
        auto d3 = d * d * d;
        std::cout << d3 << std::endl;
        real_t Nm = numberDensity() * d3() / numberDensityDPD();
        std::cout << "Coarse graining parameter: " << Nm << std::endl;
        real_t Nm_Kumar = 1.0e+8;                       // Coarse-graining parameter from Kumar2009.
        length_t dW = std::pow(Nm_Kumar * numberDensityDPD() / numberDensity(), 1.0/3.0);
        std::cout << "Characteristic length value, Nm=1.0e+08 (m):  " << dW() << std::endl;

        return Nm;
    }

    void waterInMicrochannel(const units::dpd<real_t> &dpd, const real_t& Nm) {
        density_t density{0.9970479};  // g/cm^3 at 298.15 K, Wikipedia
        density = density() * (g_to_kg / (cm_to_m * cm_to_m * cm_to_m));
        number_density_t numberDensity = 3.33324571e+28;
        //real_t Nm = 2.39993691e+12;
        real_t waterViscosity = 0.8937e-03;     // Dynamic viscosity at 298.15 K, Units are m Pa s = kg m^-1 s^-2 = N m^2 (Wikipedia)
        real_t Re = 0.33;                       // From Serhatlioglu2020, p 6936, left column. This is laminar flow.
        pressure_t pDifference = 30.0e-03;      // in bar, From Serhatlioglu2020, p 6936, left column.
        pDifference = pDifference() * 1.0e+05;  // In J/m^3 = m kg s^-2

        // Actual channel dimensions, in m.
        real_t boxX = 60.0e-06;
        real_t boxY = 50.0e-06;
        real_t boxZ = 50000.0e-06;

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "Microchannel (no units = DPD, otherwise indicated):" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        real_t r_c = dpd.length()() * nm_to_m;         // in m.
        box_ptr_t box = factory::box(boxX, boxY, boxZ);
        std::cout << "Box size (m):                                     " << *box << std::endl;
        std::cout << "Box volume (m^3):                                   " << box->volume() << std::endl;
        box = factory::box(box->lengthX() / r_c, box->lengthY() / r_c, box->lengthZ() / r_c);
        std::cout << "Box size (microchannel dimensions):               " << *box << std::endl;
        std::cout << "Box volume:                                         " << box->volume() << std::endl;
        std::cout << "Number of particles in box if number density=3.0:   " << real_t(int(3.0 * box->volume())) << std::endl;

        int N = 50000;               // Desired number of particles.
        number_density_t nd = 3.0;   // Common in DPD.
        std::cout << "DESIRED TOTAL number of DPD particles:              " << N << std::endl;
        std::cout << "DESIRED number density:                             " << nd() << std::endl;
        real_t Nl = numberDensity() * boxX * boxY / Nm;
        length_t lz = real_t(N) / Nl;
        std::cout << "Nl:                                                 " << Nl << std::endl;
        std::cout << "lz (m):                                             " << lz() << std::endl;
        box = factory::box(boxX / r_c, boxY / r_c, lz / r_c);
        std::cout << "Box dimensions:                                   " << *box << std::endl;
        std::cout << "box volume:                                         " << box->volume() << std::endl;
        std::cout << "number density:                                     " << real_t(N)/box->volume() << std::endl;
        std::cout << "Cutoff distance (m):                                " << r_c << std::endl;
        std::cout << "Cutoff distance:                                    " << 1.0 << std::endl;
        std::cout << std::endl;
        auto Dh = 2.0 * boxX * boxY / (boxX + boxY);              // Wikipedia: https://en.wikipedia.org/wiki/Hydraulic_diameter (rectangular duct)
        auto flowSpeed = Re * waterViscosity / (density() * Dh);  // Wikipedia: https://en.wikipedia.org/wiki/Reynolds_number
        std::cout << "Density (kg/m^3):                                   " << density() << std::endl;
        std::cout << "Viscosity (kg m^-1 s^-2):                           " << waterViscosity << std::endl;
        std::cout << "Dh (hydraulic parameter,L):                         " << Dh << std::endl;
        std::cout << "Reynolds number:                                    " << Re << std::endl;
        std::cout << "Fluid velocity, u (m/s):                            " << flowSpeed << std::endl;
        auto dt =  50000e-06 / flowSpeed;
        std::cout << "Time to cover 50,000 micrometers (s)                " << dt << std::endl;
        auto tau = r_c / flowSpeed;
        std::cout << "Characteristic time from flow speed (s):            " << tau << std::endl;
        auto vRMS = std::sqrt(3.0 * units::si<real_t>::kT/(dpd.mass()() * units::si<real_t>::MU));
        std::cout << "Average root mean speed (equipartition  theorem):   " << vRMS << std::endl;
        dt =  50000e-06 / vRMS;
        std::cout << "Time to cover 50,000 micrometers from RMS speed(s): " << dt << std::endl;
        real_t tauRMS = r_c / vRMS;
        std::cout << "Characteristic time from RMS speed (s):             " << tauRMS << std::endl;
        std::cout << std::endl;

        std::cout << "Pressure difference (J/m^3 = m kg s^-2):            " << pDifference() << std::endl;
        force_t F{0.0, 0.0, boxX * boxY * pDifference()};
        std::cout << "Force due to pressure gradient:                     " << norm<real_t>(F) << std::endl;
        volume_t channelVolume = boxX * boxY * boxZ;
        real_t numberOfWaterMoleculesInChannel = numberDensity() * channelVolume();
        std::cout << "Number of water molecules in channel:               " << numberOfWaterMoleculesInChannel << std::endl;
        force_t fWater{0.0, 0.0, F[2] / numberOfWaterMoleculesInChannel};
        std::cout << "Force per water molecule in channel (N):            " << fWater << std::endl;
        force_t fBead = Nm * fWater;
        std::cout << "Force per water bead (N):                           " << fBead << std::endl;
        real_t factor = tau * tau / (dpd.length()() * nm_to_m * dpd.mass()() * units::si<real_t>::MU);
        fBead *= factor;
        std::cout << "Force per water bead:                               " << fBead << std::endl;
     }
}

/**
 * Returns characteristic values for length, mass, energy, time, charge.
 * @param diameter Diameter of single water particle (in m).
 * @param temperature Temperature (in K).
 * @param rho Density at given temperature (in kg m^-3).
 * @param xy Box dimensions in x- and y-direction, respectively.
 * @param numberOfWaterFluidElements Requested number of water fluid elements.
 * @param timeStep Time step size in DPD time units,
 * @param tspan Requested time a single must cover (in s).
 */
std::tuple<real_t, real_t, real_t, real_t, real_t>
polarizableWater(const dist_t& diameter,
                 int numberOfWaterFluidElements = 6000,
                 const temperature_t& temperature = 298.15,
                 stime_t timeStep = 0.01,
                 const density_t& rho = 0.99705e+03,
                 const std::pair<real_t, real_t> xy = {6.00e-05, 5.0e-05},
                 stime_t tspan = 15) {
    util::Logger logger("simploce::polarizableWater()");
    logger.trace("Entering.");

    std::cout.setf(std::ios::scientific);

    std::cout << std::endl << std::endl;
    std::cout << "Mesoscopic Polarizable water: " << std::endl;
    std::cout << "==================================================================================="
              << "===================================" << std::endl;
    std::cout << "A single fluid water element is represented by TWO (2) connected beads (i.e. particles): " << std::endl;
    std::cout << "(1) Heavier central water (CW)," << std::endl;
    std::cout << "(2) Lighter dipolar (DP) bead." << std::endl;
    std::cout << std::endl;

    const mass_t massWaterMolecule{2.99063e-26};             // kg.

    const number_density_t ndDPDStandard{3.0};               // Standard value for the TOTAL bead number
                                                             // density (Groot and Warren, 1997), in DPD units.

    const number_density_t ndDPD = 0.5 * ndDPDStandard();    // TOTAL number of beads (i.e. particles) / volume. The
                                                             // factor 0.5 is because each fluid element is
                                                             // represented by TWO beads instead of just one.
                                                             // Thus, this number refers to each of the two TYPES of
                                                             // beads, that is, the CW and DP bead.

    dist_t lx{xy.first};                                     // Box size in x-direction, in m.
    dist_t ly{xy.second};                                    // Box size in y-direction, in m.

    // Requested values/setup.
    std::cout << "Requested setup:" << std::endl;
    std::cout << diameter << ": Diameter water bead (in m), -if- a single water fluid element is represented "
                             "by a single DPD bead." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << numberOfWaterFluidElements;
    std::cout << ": Number of water fluid elements." << std::endl;
    std::cout << temperature << ": Temperature (in K)." << std::endl;
    std::cout << rho << ": Water density (in kg m^-3) at given temperature." << std::endl;
    std::cout << lx << ": Box dimension in x-direction (in m)." << std::endl;
    std::cout << ly << ": Box dimension in y-direction (in m)." << std::endl;
    std::cout << tspan << ": Required time span attainable (in s)." << std::endl;

    int N = 2 * numberOfWaterFluidElements; // TOTAL number of beads. Each fluid element is represented by 2 beads.
    std::cout << std::setw(conf::REAL_WIDTH) << N;
    std::cout << ": TOTAL number of beads (= 2 x Number of water fluid elements)." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << N/2 << ": Number of CW beads." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << N/2 << ": Number of DP beads." << std::endl;

    dist_t dW{diameter};                                                // In m.
    temperature_t Tref{temperature};                                    // in K.
    energy_t kB_Tref = units::si<real_t>::KB * Tref();                  // in J.
    number_density_t ndWaterMolecules = rho() / massWaterMolecule();    // Number density water molecules, in m^-3.
    real_t Nm =  2.0 * ndWaterMolecules() * dW() * dW() * dW() / ndDPDStandard(); // Coarse-graining parameter.
    mass_t mW = massWaterMolecule * Nm;
    std::cout << ndDPDStandard << ": Standard DPD particle number density." << std::endl;
    std::cout << massWaterMolecule << ": Mass single water molecule (in kg)." << std::endl;
    std::cout << ndWaterMolecules << ": Number density (in m^-3) of water molecules." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << Nm << ": Coarse-graining parameter." << std::endl;
    std::cout << mW << ": Mass of a single water fluid element." << std::endl;
    std::cout << std::endl;

    // Box dimension in z-direction.
    std::cout << "Box dimension in z-direction:" << std::endl;
    real_t volume_lz = lx() * ly();                                  // Volume * lz, in m^2.
    real_t numberWaterMolecules_lz = ndWaterMolecules() * volume_lz; // (Number of water molecuĺes) * lz.
    real_t M_lz = numberWaterMolecules_lz / Nm;                      // (Number of water fluid elements) * lz.
    real_t N_lz = 2.0 * M_lz;                                        // (Total number of beads) * lz.
    dist_t lz = N / N_lz;                   // Box dimension in z-direction.
    std::cout << std::setw(conf::REAL_WIDTH) << volume_lz << ": Volume * lz (in m^2)." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << numberWaterMolecules_lz << ": (Number of water molecules) * lz" << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << M_lz << ": (Number of water fluid elements) * lz" << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << N_lz << ": (Total number of beads) * lz" << std::endl;
    std::cout << lz << ": Box dimension (lz) in z-direction (in m)" << std::endl;
    std::cout << lz/1.0e-06 << ": Box dimension (lz) in z-direction (in micrometer)" << std::endl;
    std::cout << dist_t{lz() / dW()} << ": Box dimension (lz) in z-direction (in DPD units)" << std::endl;
    std::cout << lx() / dW << ": Box dimension (lx) in the x-direction." << std::endl;
    std::cout << ly() / dW << ": Box dimension (ly) in the y-direction." << std::endl;
    volume_t volumeDPD = lx() * ly() * lz() / (dW() * dW() * dW());        // Volume in DPD units.
    std::cout << volumeDPD << ": Volume (in DPD units)." << std::endl;
    number_density_t ndActual{N/volumeDPD()};
    std::cout << ndActual << ": Total bead number density (DPD units), -SHOULD BE 3.0-." << std::endl;
    std::cout << std::endl;

    // Time.
    std::cout << "Time:" << std::endl;
    dist_t cutoff = dW;
    mass_t mass = mW;
    real_t vref = std::sqrt(kB_Tref() / mass()); // Thermal velocity in one dimension.
    stime_t tauW = cutoff() / vref;
    stime_t dt_in_s = timeStep() * tauW();
    std::cout << cutoff << ": Cutoff distance (in m)." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << vref << ": Thermal speed in one dimension (in m s^-1)." << std::endl;
    std::cout << timeStep << ": Time step (in DPD units)." << std::endl;
    std::cout << dt_in_s << ": Time step (in s)." << std::endl;
    int n = int(tspan()/dt_in_s());
    std::cout << std::setw(conf::REAL_WIDTH) << n << ": Number of steps to cover time span of " << tspan() << " (in s)." << std::endl;
    n = int((lz() / vref) / dt_in_s());
    std::cout << std::setw(conf::REAL_WIDTH) << n << ": Number of steps to cover ";
    std::cout << std::setprecision(2) << lz()/1.0e-06 << " (lz, in micrometer)." << std::endl;
    std::cout << std::endl;

    // Charge
    std::cout << "Charge:" << std::endl;
    charge_t qW = std::sqrt(units::si<real_t>::E0 * dW() * kB_Tref());  // In C
    charge_t one_C_in_dpd_units = 1.0 * units::si<real_t>::E / qW();
    charge_t charge_CW_Riniker_dpd_units = 0.575 * units::si<real_t>::E / qW();
    std::cout << one_C_in_dpd_units << ": 1 Coulomb in DPD units." << std::endl;
    std::cout << charge_CW_Riniker_dpd_units << ": Charge value of CW/DP (0.575 e) in Riniker2011 in DPD units." << std::endl;
    std::cout << std::endl;

    // Halve attractive quartic bonded potential, U(r) = 0.5 * k * (r - r0)^4
    std::cout << "Halve attractive quartic bonded potential, U(r) = 0.5 * k * (r - r0)^4" << std::endl;
    real_t kq = 2.0e+06; // force constant, in kJ mol^-1 nm^-4.
    std::cout << std::setw(conf::REAL_WIDTH) << kq << ": Force constant in kJ mol^-1 nm^-4." << std::endl;
    kq *= 1000.0;                                      // in J mol^-1 nm^-4
    kq /= (nm_to_m * nm_to_m * nm_to_m * nm_to_m);     // in J mol^-1 m^-4
    kq /= units::si<real_t>::NA;            // in J m^-4
    std::cout << std::setw(conf::REAL_WIDTH) << kq << ": Force constant in J m^-4." << std::endl;
    kq *= kB_Tref() / (dW() * dW() * dW() * dW());
    std::cout << std::setw(conf::REAL_WIDTH) << kq << ": Force constant in DPD units." << std::endl;
    std::cout << std::endl;

    // Calculating charge value for CW using halve attractive harmonic potential, U(r) = 0.5 * (r-r0)^2.
    // Already in DPD units.
    std::cout << "Halve attractive quartic harmonic potential, U(r) = 0.5 * k * (r - r0)^2" << std::endl;
    auto eps_r_W = 78.0;
    auto alpha = 0.5 * math::constants<real_t>::PI * (eps_r_W - 1.0) / (eps_r_W + 2);
    std::cout << std::setw(conf::REAL_WIDTH) << eps_r_W << ": Dielectric constant of water." << std::endl;
    std::cout << std::setw(conf::REAL_WIDTH) << alpha << ": Polarizability of sphere of radius=0.5 l." << std::endl;
    std::cout << "Force constant HP, Charge CW (in DPD units)" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    for (int k = 1; k != 10; ++k) {
        real_t fc = k * 100.0;
        charge_t Q_CW = std::sqrt(alpha * fc);
        std::cout << std::setw(conf::REAL_WIDTH) << fc << " " << std::setw(conf::REAL_WIDTH) << Q_CW << std::endl;
    }
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << std::endl;

    // Display characteristic values in SI units.
    std::cout << std::endl;
    std::cout << "Characteristic values in SI:" << std::endl;
    std::cout << dW << ": Characteristic value for length (in m)." << std::endl;
    std::cout << kB_Tref << ": Characteristic value for energy (in J)." << std::endl;
    std::cout << mW << ": Characteristic value for mass (in kg)." << std::endl;
    std::cout << tauW << ": Characteristic value for time (in s)." << std::endl;
    std::cout << qW << ": Characteristic value for charge (in C)." << std::endl;
    std::cout << "==================================================================================="
              << "===================================" << std::endl;
    std::cout << std::endl << std::endl;

    logger.trace("Leaving.");
    return std::make_tuple(dW(), mW(), kB_Tref(), tauW(), qW());
}

int main(int argc, char *argv[]) {
    util::Logger::changeLogLevel(util::Logger::LOGINFO);
    util::Logger logger("simploce::s-dpd-units::main");

    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    real_t diameter{6.0e-06};        // m.
    real_t temperature{298.15};      // K.
    real_t rhoWater{0.99705e+03};    // kg m^-3
    int numberOfWaterFluidElements{6000};
    real_t timeStep{0.01};           // in DPD units.

    po::options_description usage("Usage");
    usage.add_options() (
        "fn-particle-spec-catalog,s",
        po::value<std::string>(&fnParticleSpecCatalog),
        "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
    )(
        "temperature,T",
        po::value<real_t>(&temperature),
        "Temperature (K). Default is 298.15 K."
    )(
        "number-of-water-fluid-elements,M",
        po::value<int>(&numberOfWaterFluidElements),
        "Number of water fluid elements. Default is 6000."
    )(
        "timestep",
        po::value<real_t>(&timeStep),
        "Time step in DPD units. Default value is 0.01."
    )
    (
        "verbose,v",
        "Verbose"
    )(
         "help,h",
         "Help message"
    );

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, usage), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << "Compute characteristic values for DPD" << std::endl;
        std::cout << usage << "\n";
        return 0;
    }
    if (vm.count("fn-particle-spec-catalog,s")) {
        fnParticleSpecCatalog = vm["fn-particle-spec-catalog"].as<std::string>();
    }
    if (vm.count("temperature")) {
        temperature = vm["temperature"].as<real_t>();
    }
    if (vm.count("number-of-water-fluid-elements")) {
        numberOfWaterFluidElements = vm["number-of-water-fluid-elements"].as<int>();
    }
    if (vm.count("timestep")) {
        timeStep = vm["timestep"].as<real_t>();
    }
    if (vm.count("verbose") ) {
        util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    }

    std::cout.setf(std::ios::scientific);
    std::cout.precision(5);
    auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
    //auto dpd = polystyreneInWaterChannel(catalog);
    //auto Nm = water(dpd, catalog);
    //waterInMicrochannel(dpd, Nm);

    polarizableWater(diameter, numberOfWaterFluidElements, temperature, timeStep);
}
