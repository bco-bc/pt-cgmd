/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 12:50 PM
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/p-util.hpp"
#include "simploce/particle/atomic-content-handler.hpp"
#include "simploce/chem/input-source.hpp"
#include "simploce/chem/pdb-reader.hpp"
#include "simploce/util/util.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include <stdexcept>
#include <cmath>
#include <utility>
#include <random>

namespace simploce {

    ParticleSystemFactory::ParticleSystemFactory(spec_catalog_ptr_t catalog) :
        catalog_{std::move(catalog)}
    {
        if ( !catalog_ ) {
            util::Logger logger("simploce::ParticleSystemFactory::ParticleSystemFactory");
            util::logAndThrow(logger,
                              "simploce::ParticleSystemFactory: Missing particle specifications catalog.");
        }
    }

    p_system_ptr_t
    ParticleSystemFactory::diatomic(const dist_t& distance,
                                    const spec_ptr_t& spec,
                                    const temperature_t& temperature) {
        util::Logger logger{"simploce::ParticleSystemFactory::diatomic"};

        // Atomistic model.
        auto atomistic = factory::atomistic();
        std::vector<p_ptr_t> atoms{};

        // Atom #1
        position_t r1{-0.5 * distance() - 0.03, 0.0, 0.0};
        std::string name1 = spec->name() + "1";
        auto atom_1 = atomistic->addAtom(name1, spec);
        atom_1->position(r1);
        util::assignVelocity(atom_1, temperature);

        // Motion along the x-axis.
        auto v = atom_1->velocity();
        v[1] = 0.0;
        v[2] = 0.0;
        atom_1->velocity(v);

        atoms.emplace_back(atom_1);

        // Atom #2
        position_t r2{0.5 * distance() + 0.03, 0.0, 0.0};
        std::string name2 = spec->name() + "2";
        auto atom_2 = atomistic->addAtom(name2, spec);
        atom_2->position(r2);
        util::assignVelocity(atom_2, temperature);

        // Initially opposite motion.
        v = -1.0 * atom_1->velocity();
        atom_2->velocity(v);

        atoms.emplace_back(atom_2);

        // Make bond.
        id_pair_t bond = std::make_pair(atom_1->id(), atom_2->id());
        std::vector<id_pair_t> bonds{bond};
        auto atomGroup = ParticleGroup::make(atoms, bonds);
        atomistic->addParticleGroup(atomGroup);

        // Place diatomic in a very large box.
        box_ptr_t box(new Cube<real_t>{1000.0});
        atomistic->box(box);

        util::removeOverallLinearMomentum(atomistic);

        // Done.
        return std::move(atomistic);
    }

    p_system_ptr_t
    ParticleSystemFactory::polarizableWater(const box_ptr_t& box,
                                            std::size_t nLimit,
                                            const density_t& densitySI,
                                            const temperature_t& temperature) {
        const dist_t DISTANCE_CW_DP = 0.2; // nm.

        util::Logger logger{"ParticleSystemFactory::polarizableWater"};

        logger.info("Creating coarse grained particle model for polarizable water.");
        logger.info("Temperature: " + boost::lexical_cast<std::string>(temperature));

        // Convert kg/m^3 to u/nm^3.
        density_t density = densitySI / (units::si<real_t>::MU * 1.0e+27);
        logger.info("Requested density (kg/m^3): " + boost::lexical_cast<std::string>(densitySI));
        logger.info("Requested density (u/nm^3): " + boost::lexical_cast<std::string>(density));

        // Box details.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        logger.info("Box size (nm): " + boost::lexical_cast<std::string>(box->size()));
        logger.debug("Box volume (nm^3): " + boost::lexical_cast<std::string>(volume));

        // How many CG particles? A single polarizable CG water particle represents
        // 5 water molecules. Each CG water particles consists of two connected
        // CG particles (CW and DP).
        auto h2oSpec = catalog_->molecularWater();
        auto cwSpec = catalog_->lookup("CW");
        auto dpSpec = catalog_->lookup("DP");
        logger.debug(
                "\"Ideal\" distance between CW and DP (nm): " + boost::lexical_cast<std::string>(DISTANCE_CW_DP)
        );

        number_density_t atNumberDensity = density / h2oSpec->mass();
        int natWaters = int(atNumberDensity * volume);
        int ncgWaters = int(real_t(natWaters) / 5.0);
        if ( ncgWaters > int(nLimit) ) {
            ncgWaters = int(nLimit);
            logger.warn("Number of requested CG waters is higher than the maximum number of CG water allowed.");
        }
        logger.info("Generating " + util::toString(ncgWaters) + " CG water.");
        number_density_t cgNumberDensity = real_t(ncgWaters) / volume();
        std::size_t nOneDirection = int(std::pow(ncgWaters, 1.0/3.0));
        density_t cgDensity = cgNumberDensity * (cwSpec->mass() + dpSpec->mass());
        logger.debug("CG: Requested density (u/nm^3): " +  boost::lexical_cast<std::string>(cgDensity));
        logger.debug("AT: Requested number density (1/nm^3): " + boost::lexical_cast<std::string>(atNumberDensity));
        logger.debug("AT: Requested number density (1/m^3): " + boost::lexical_cast<std::string>(atNumberDensity*1.0e+27));
        logger.debug("CG: Requested number density (1/nm^3): " + boost::lexical_cast<std::string>(cgNumberDensity));
        logger.debug("CG: Requested number density (1/m^3): " + boost::lexical_cast<std::string>(cgNumberDensity*1.0e+27));
        logger.debug("AT: Required number of water molecules: " + boost::lexical_cast<std::string>(natWaters));
        logger.info("CG: Requested number of water water groups: " + boost::lexical_cast<std::string>(ncgWaters));

        // Spacing between DW water beads.
        const length_t spacing{0.53};  // Roughly the location of the first peak of g(r),
                                          // in nm, Figure 5 in Riniker et al.
        logger.debug("Spacing between CW beads: " + boost::lexical_cast<std::string>(spacing));

        std::size_t nx = util::nint(Lx / spacing);
        std::size_t ny = util::nint(Ly / spacing);
        std::size_t nz = util::nint(Lz / spacing);

        logger.debug("Number of coordinates in x-direction: " + boost::lexical_cast<std::string>(nx));
        logger.debug("Number of coordinates in y-direction: " + boost::lexical_cast<std::string>(ny));
        logger.debug("Number of coordinates in z-direction: " + boost::lexical_cast<std::string>(nz));
        logger.debug("Number of coordinates in any direction: " + boost::lexical_cast<std::string>(nOneDirection));
        logger.debug("Can create a maximum of " + util::toString(nx * ny * nz) + " CW beads.");

        // Start from an empty particle model.
        auto coarseGrained = factory::coarseGrained();

        // Add beads.
        real_t x0 = 0.01; //-0.5 * Lx();
        real_t y0 = 0.01; //-0.5 * Ly();
        real_t z0 = 0.01; //-0.5 * Lz();
        std::size_t counter = 0;
        int i = 0, j = 0, k = 0;
        while ( i < nx && counter < ncgWaters ) {
            //real_t x = x0 + (i + 0.5) * spacing();
            real_t x = x0 + i * spacing();
            while ( j < ny && counter < ncgWaters ) {
                //real_t y = y0 + (j + 0.5) * spacing();
                real_t y = y0 + j * spacing();
                while ( k < nz && counter < ncgWaters ) {
                    //real_t z = z0 + (k + 0.5) * spacing();
                    real_t z = z0 + k * spacing();

                    // Single water group.
                    std::vector<p_ptr_t> beads{};
                    std::vector<id_pair_t> bonds{};

                    position_t r1{x,y,z};
                    std::string index =  boost::lexical_cast<std::string>(counter + 1);
                    std::string name = "CW" + index;
                    auto cwBead = coarseGrained->addBead(name, cwSpec);
                    cwBead->position(r1);
                    util::assignVelocity(cwBead, temperature);
                    beads.emplace_back(cwBead);

                    // Place the DP parallel to one of the 3 coordinate axes.
                    position_t r2{};
                    int l = int(util::random() * 3.0);
                    switch (l) {
                        case 0: {
                            real_t coord = x + DISTANCE_CW_DP();
                            r2 = position_t{coord, y, z};
                            break;
                        }
                        case 1: {
                            real_t coord = y + DISTANCE_CW_DP();
                            r2 = position_t{x,coord,z};
                            break;
                        }
                        default: {
                            real_t coord = z + DISTANCE_CW_DP();
                            r2 = position_t{x,y,coord};
                            break;
                        }
                    }
                    name = "DP" + index;
                    auto dpBead = coarseGrained->addBead(name, dpSpec);
                    dpBead->position(r2);
                    util::assignVelocity(dpBead, temperature);
                    beads.emplace_back(dpBead);

                    id_pair_t pair = std::make_pair(cwBead->id(), dpBead->id());
                    bonds.emplace_back(pair);

                    coarseGrained->addParticleGroup(beads, bonds);

                    counter += 1;
                    k += 1;
                }
                k = 0;
                j += 1;
            }
            j = 0;
            i += 1;
        }
        logger.info("Created " + util::toString(coarseGrained->numberOfBeads()) + " beads.");
        logger.info("Created " + util::toString(coarseGrained->numberOfParticleGroups()) + " water groups.");
        if ( coarseGrained->numberOfParticleGroups() < ncgWaters ) {
            logger.warn("The number of created water groups is less then the number of requested water groups.");
        }
        coarseGrained->box(box);

        util::removeOverallLinearMomentum(coarseGrained);

        return std::move(coarseGrained);
    }

    p_system_ptr_t
    ParticleSystemFactory::simpleElectrolyte(const box_ptr_t& box,
                                             const molarity_t& molarity,
                                             const temperature_t& temperature)
    {
        util::Logger logger{"ParticleSystemFactory::simpleElectrolyte"};

        logger.info("Creating atomistic particle model for a simple electrolyte solution:");

        logger.info("Molarity (mol/l): " + boost::lexical_cast<std::string>(molarity()));
        logger.info("Temperature (K) : " + boost::lexical_cast<std::string>(temperature()));

        spec_ptr_t NA = catalog_->lookup("Na+");
        spec_ptr_t CL = catalog_->lookup("Cl-");

        // Box details.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        logger.info("Box size (nm): " + boost::lexical_cast<std::string>(box->size()));
        logger.debug("Box volume (nm^3): " + boost::lexical_cast<std::string>(volume));

        // Determine required number of ions.
        number_density_t numberDensity =
            molarity() * units::si<real_t>::NA / units::mu<real_t>::l_to_nm3;
        int nIons = util::nint(2.0 * numberDensity * volume);
        nIons = (nIons % 2 != 0 ? nIons + 1 : nIons);  // Neutrality condition.
        logger.info("TOTAL number of ions requested: " + boost::lexical_cast<std::string>(nIons));

        // Spacing between particles.
        const length_t spacing = 0.6; // The location of the first peak of g(r) for
                                      // 0.1 M NaCl is around 0.3 nm.
        std::size_t nx = util::nint(Lx / spacing);
        std::size_t ny = util::nint(Ly / spacing);
        std::size_t nz = util::nint(Lz / spacing);
        logger.debug("Distance spacing between ion particles: " + boost::lexical_cast<std::string>(spacing));
        logger.debug("Number of coordinates in x-direction: " +  boost::lexical_cast<std::string>(nx));
        logger.debug("Number of coordinates in y-direction: " +  boost::lexical_cast<std::string>(ny));
        logger.debug("Number of coordinates in z-direction: " +  boost::lexical_cast<std::string>(nz));

        // Start from an empty particle model.
       auto atomistic = factory::atomistic();

        // Create ions.
        std::size_t counter = 0;
        int i = 0, j = 0, k = 0;

        while ( i < nx && counter < nIons ) {
            while ( j < ny && counter < nIons ) {
                while ( k < nz && counter < nIons ) {

                    // Position.
                    real_t x = (i + 0.5) * spacing() + util::random() * 0.1;
                    real_t y = (j + 0.5) * spacing() - util::random() * 0.1;
                    real_t z = (k + 0.5) * spacing() + util::random() * 0.1;
                    position_t r{x,y,z};
                    if ( counter % 2 == 0 ) {
                        // Na+
                        std::string name = NA->name() + util::toString(counter + 1 );
                        auto ion = atomistic->addAtom(name, NA);
                        ion->position(r);
                        util::assignVelocity(ion, temperature);
                    } else {
                        // Cl-
                        std::string name = CL->name() + util::toString(counter + 1);
                        auto ion = atomistic->addAtom(name, CL);
                        ion->position(r);
                        util::assignVelocity(ion, temperature);
                    }

                    counter += 1;
                    k += 1;
                }
                k = 0;
                j += 1;
            }
            j = 0;
            i += 1;
        }
        atomistic->box(box);

        logger.info("Number of ion particles generated: " + boost::lexical_cast<std::string>(counter));
        if (counter < nIons) {
            logger.warn("Number of generated ions is smaller than requested.");
        }
        logger.info("Ion number density (1/nm^3): " + boost::lexical_cast<std::string>(real_t(counter)/volume()));

        util::removeOverallLinearMomentum(atomistic);

        return std::move(atomistic);
    }

    p_system_ptr_t
    ParticleSystemFactory::argon(const box_ptr_t& box,
                                 std::size_t nLimit,
                                 const density_t& densitySI,
                                 const temperature_t& temperature)
    {
        util::Logger logger{"simploce::ParticleSystemFactory::argon"};

        logger.info("Creating atomistic particle model for argon.");
        logger.info("Requested temperature: " + util::toString(temperature()));

        // Conversion of kg/m^3 to u/nm^3.
        density_t density = densitySI / (units::si<real_t>::MU * 1.0e+27);

        // Argon specification.
        auto spec = catalog_->lookup("Ar");
        logger.info("Argon specification:");
        util::log(logger, *spec);

        number_density_t numberDensity = density / spec->mass();
        logger.info("Requested density (kg/m^3): " + util::toString(densitySI));
        logger.info("Requested density (u/nm^3): " + util::toString(density));
        logger.info("Requested number density (1/nm^3):" + util::toString(numberDensity));

        // Spacing between argon atoms.
        const length_t spacing{0.3};
        logger.debug("Spacing between LJ atoms: " + util::toString(spacing));

        // Box dimensions.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        logger.info("Requested box size (nm): " + util::toString(box->size()));
        logger.debug("Box volume (nm^3): " + util::toString(volume));

        int nAtoms = int(numberDensity() * volume());
        logger.info("Requested number of argon atoms: " + util::toString(nAtoms));
        if (nAtoms > int(nLimit)) {
            logger.warn("Number of requested atoms is higher that maximum number of atoms allowed.");
            nAtoms = int(nLimit);
        }
        logger.info("Generating " + util::toString(nAtoms) + " argon atoms.");

        int nx = util::nint(Lx / spacing);
        int ny = util::nint(Ly / spacing);
        int nz = util::nint(Lz / spacing);
        logger.debug("Number of coordinates in x-direction: " + util::toString(nx));
        logger.debug("Number of coordinates in y-direction: " + util::toString(ny));
        logger.debug("Number of coordinates in z-direction: " + util::toString(nz));

        // Start from an empty particle model.
        auto atomistic = factory::atomistic();

        // Add atoms.
        real_t x0 = 0.0;
        real_t y0 = 0.0;
        real_t z0 = 0.0;
        int counter = 0;
        int i = 0, j = 0, k = 0;
        while ( i < nx && counter < nAtoms ) {
            real_t x = x0 + (i + 0.5) * spacing();
            while ( j < ny && counter < nAtoms ) {
                real_t y = y0 + (j + 0.5) * spacing();
                while ( k < nz && counter < nAtoms ) {
                    real_t z = z0 + (k + 0.5) * spacing();
                    position_t r{x,y,z};
                    std::string name = spec->name() + util::toString(counter + 1);
                    auto atom = atomistic->addAtom(name, spec);
                    atom->position(r);
                    util::assignVelocity(atom, temperature);
                    counter += 1;
                    k += 1;
                }
                k = 0;
                j += 1;
            }
            j = 0;
            i += 1;
        }
        atomistic->box(box);
        logger.info("Number of Argon atoms generated: " + util::toString(atomistic->numberOfAtoms()));
        numberDensity = real_t(atomistic->numberOfAtoms()) / volume();
        logger.info("Number density (1/nm^3): " + util::toString(numberDensity));
        density = spec->mass()() * numberDensity();
        logger.info("Density (u/nm3): " + util::toString(density));
        density_t densSI = density() * units::si<real_t>::MU * 1.0e+27;
        logger.info("Density (kg/m^3): " + util::toString(densSI));

        util::removeOverallLinearMomentum(atomistic);

        return std::move(atomistic);
    }

    p_system_ptr_t
    ParticleSystemFactory::polymerSolution(const box_ptr_t& box,
                                           int chainLength,
                                           const std::string& monomericUnitBeadSpecName,
                                           int numberOfPolymers,
                                           const length_t& spacing,
                                           int numberOfWaters,
                                           const std::string& waterBeadSpecName,
                                           const temperature_t& temperature,
                                           bool placeRandom) {
        util::Logger logger("simploce::ParticleSystemFactory::polymerSolution");
        logger.trace("Entering.");

        logger.info("Creating coarse-grained particle system of a polymer solution.");
        logger.info(util::toString(*box) + ": Requested box dimensions.");
        logger.info(std::to_string(chainLength) + ": Requested polymer chain length (number of monomeric units).");
        logger.info(monomericUnitBeadSpecName + ": Monomeric unit particle specification name.");
        logger.info(std::to_string(numberOfPolymers) + ": Requested number of polymers.");
        logger.info(std::to_string(spacing()) + ": Requested spacing between monomeric units.");
        logger.info(std::to_string(numberOfWaters) + ": Requested number of water beads.");
        logger.info(waterBeadSpecName + ": Water specification name.");
        logger.info(std::to_string(temperature()) + ": Requested temperature.");
        auto numberOfBeads = numberOfWaters + numberOfPolymers * chainLength;
        logger.info(std::to_string(numberOfPolymers * chainLength) + ": Number of polymer beads required.");
        logger.info(std::to_string(numberOfWaters) + ": Number of water beads required.");
        logger.info(std::to_string(numberOfBeads) + ": Total number of beads required.");
        logger.info(std::to_string(placeRandom) + ": Place beads randomly in box.");

        auto spec = catalog_->lookup(monomericUnitBeadSpecName);
        auto waterSpec = catalog_->lookup(waterBeadSpecName);
        if (box->size() < (chainLength - 1) * spacing() + 2.0 * spec->radius()() && numberOfPolymers > 0) {
            util::logAndThrow(logger, "Box size is too small for required polymer chain length.");
        }
        length_t longestSide{box->lengthX()};
        std::size_t longestSideIndex = 0;
        for (std::size_t k = 1; k != 3; ++k) {
            auto length = box->operator[](k);
            if ( box->operator[](k) > longestSide() ) {
                longestSideIndex = k;
                longestSide = length;
            }
        }
        logger.debug(std::to_string(longestSide()) + ": Longest box side.");
        logger.debug(std::to_string(longestSideIndex) + ": Longest box side index.");
        if (longestSideIndex != 2 && numberOfPolymers > 0) {
            util::logAndThrow(logger, "Longest box side should be in the z-direction.");
        }

        auto polymerSystem = factory::coarseGrained();
        auto boxX = box->lengthX();
        auto boxY = box->lengthY();
        auto boxZ = box->lengthZ();

        if (placeRandom) {
            int counter{0};

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<real_t> dis(0.0, 1.0);

            // Polymers. Placed randomly in the box.
            for (auto i = 0; i != numberOfPolymers; ++i) {
                std::vector<p_ptr_t> particles;
                auto x = dis(gen) * boxX;
                auto y = dis(gen) * boxY;
                auto z = dis(gen) * boxZ;
                position_t r{x, y, z};
                for (auto j = 0; j != chainLength; ++j) {
                    r[2] = z + real_t(j) * spacing();
                    auto name = spec->name() + util::toString(counter);
                    auto particle = polymerSystem->addBead(name, spec);
                    particle->position(r);
                    util::assignVelocity(particle, temperature, true);
                    particles.emplace_back(particle);
                    counter += 1;
                }
                std::vector<id_pair_t> bonds;
                for (std::size_t j = 0; j != (chainLength - 1); ++j) {
                    auto id1 = particles[j]->id();
                    auto id2 = particles[j + 1]->id();
                    bonds.emplace_back(id_pair_t{id1, id2});
                }
                polymerSystem->addParticleGroup(particles, bonds);
            }

            // Add water beads, placed randomly in the box.
            for (auto i = 0; i != numberOfWaters; ++i) {
                auto x = dis(gen) * boxX;
                auto y = dis(gen) * boxY;
                auto z = dis(gen) * boxZ;
                position_t r{x, y, z};
                auto name = waterBeadSpecName + util::toString(counter);
                auto particle = polymerSystem->addBead(name, waterSpec);
                particle->position(r);
                util::assignVelocity(particle, temperature, true);
                counter += 1;
            }
            logger.info(std::to_string(counter) + ": Number of beads placed randomly in box.");
        } else {
            // Initialize counters.
            int counter{0};
            int maxCounterX = int(box->lengthX() / spacing());
            int maxCounterY = int(box->lengthY() / spacing());
            int maxCounterZ = int(box->lengthZ() / spacing());
            logger.debug(std::to_string(maxCounterX) + ": Maximum value counterX.");
            logger.debug(std::to_string(maxCounterY) + ": Maximum value counterY.");
            logger.debug(std::to_string(maxCounterZ) + ": Maximum value counterZ.");
            auto numberOfGridPositions = maxCounterX * maxCounterY * maxCounterZ;
            logger.debug(std::to_string(numberOfGridPositions) + ": Number of available grid positions.");
            int counterX = 0;
            int counterY = 0;
            int counterZ = 0;

            for (std::size_t i = 0; i != numberOfPolymers; ++i) {
                // Create a polymer.
                std::vector<p_ptr_t> particles;
                for (std::size_t j = 0; j != chainLength; ++j) {
                    //std::clog << "Polymer: " << counter << " " << counterX << " " << counterY << " " << counterZ << std::endl;
                    auto x = counterX * spacing();
                    auto y = counterY * spacing();
                    auto z = counterZ * spacing();
                    position_t r{x, y, z};
                    auto name = spec->name() + util::toString(counter);
                    auto particle = polymerSystem->addBead(name, spec);
                    counter += 1;
                    particle->position(r);
                    util::assignVelocity(particle, temperature, true);
                    particles.emplace_back(particle);
                    counterZ += 1;
                    if (counterZ == maxCounterZ) {
                        logger.debug(std::to_string(i) + ": Number of polymers added.");
                        util::logAndThrow(logger,
                                          "Box length in z-direction is too short to accommodate "
                                          "more polymer chains. Try again with reduced spacing.");
                    }
                }
                std::vector<id_pair_t> bonds;
                for (std::size_t j = 0; j != (chainLength - 1); ++j) {
                    auto id1 = particles[j]->id();
                    auto id2 = particles[j + 1]->id();
                    bonds.emplace_back(id_pair_t{id1, id2});
                }
                polymerSystem->addParticleGroup(particles, bonds);

                if (counter != numberOfPolymers * chainLength) {
                    counterZ -= chainLength;
                    counterX += 1;
                    if (counterX == maxCounterX) {
                        counterX = 0;
                        counterY += 1;
                        if (counterY == maxCounterY) {
                            counterZ += chainLength;
                            counterY = 0;
                        }
                    }
                }
            }
            counterX = 0;
            counterY = 0;
            for (auto i = 0; i != numberOfWaters; ++i) {
                //std::clog << "Water: " << counter << " " << counterX << " " << counterY << " " << counterZ << std::endl;
                auto x = counterX * spacing();
                auto y = counterY * spacing();
                auto z = counterZ * spacing();
                auto name = waterBeadSpecName + util::toString(counter);
                position_t r{x, y, z};
                auto particle = polymerSystem->addBead(name, waterSpec);
                util::assignVelocity(particle, temperature, true);
                particle->position(position_t{x, y, z});
                counter += 1;

                counterX += 1;
                if (counterX == maxCounterX) {
                    counterX = 0;
                    counterY += 1;
                    if (counterY == maxCounterY) {
                        counterY = 0;
                        counterZ += 1;
                        if (counterZ == maxCounterZ) {
                            util::logAndThrow(logger,
                                              "Box too small to accommodate "
                                              "more water. Try again with reduced spacing.");
                        }
                    }
                }
            }
            logger.info(std::to_string(counter) + ": Number of beads placed on a grid in the box.");
        }
        util::removeOverallLinearMomentum(polymerSystem);

        // Box
        polymerSystem->box(box);

        logger.info(std::to_string(polymerSystem->numberOfBeads()) + ": Total number of beads created.");
        logger.info(std::to_string(polymerSystem->numberOfParticleGroups()) + ": Number of particle groups created.");
        logger.info(std::to_string(polymerSystem->numberOfFreeParticles()) + ": Number of free beads created." );

        logger.trace("Leaving");
        return std::move(polymerSystem);
    }

    void
    ParticleSystemFactory::addParticleBoundary(const p_system_ptr_t& particleSystem,
                                               dist_t spacing,
                                               Plane plane,
                                               bool excludeCorner) {
        util::Logger logger{"simploce::ParticleSystemFactory::addParticleBoundary"};
        logger.info("Adding boundary particles.");
        logger.info("Spacing (nm): " + util::toString(spacing()));
        logger.info("Plane: " + plane.toString());
        auto box = particleSystem->box();
        /*auto size = box->size();
        size += spacing();
        box = factory::box(size);
        particleSystem->box(box);*/

        length_t limit1, limit2, dis1, dis2;

        if (plane == Plane::XY) {
            int n = int(box->lengthX() / spacing());
            dis1 = box->lengthX() / n;
            n = int(box->lengthY() / spacing());
            dis2 = box->lengthY() / n;
            limit1 = box->lengthX();
            limit2 = box->lengthY();
        } else if (plane == Plane::YZ) {
            int n = int(box->lengthY() / spacing());
            dis1 = box->lengthY() / n;
            n = int(box->lengthZ() / spacing());
            dis2 = box->lengthZ() / n;
            limit1 = box->lengthY();
            limit2 = box->lengthZ();
        } else {
            int n = int(box->lengthZ() / spacing());
            dis1 = box->lengthZ() / n;
            n = int(box->lengthX() / spacing());
            dis2 = box->lengthX() / n;
            limit1 = box->lengthZ();
            limit2 = box->lengthX();
        }
        logger.info("Upper length limit #1" + util::toString(limit1));
        logger.info("Upper length limit #2" + util::toString(limit2));
        std::vector<p_ptr_t> particles;
        std::vector<id_pair_t> bonds;  // No bonds.

        real_t coord1 = 0.0;
        real_t coord2 = excludeCorner ? dis2() : 0.0;
        limit2 = excludeCorner ? limit2 -= dis2 : limit2;

        int counter = 0;
        do {
            do {
                counter += 1;
                std::string name = "BP" + util::toString(counter);
                auto p1 = particleSystem->addParticle(name,catalog_->staticBP());
                if (plane == Plane::XY) {
                    position_t r{coord1, coord2, 0.0};
                    p1->position(r);
                } else if (plane == Plane::YZ) {
                    position_t r{0.0, coord1, coord2};
                    p1->position(r);
                } else {
                    position_t r{coord2, 0.0, coord1};
                    p1->position(r);
                }
                particles.emplace_back(p1);
                counter += 1;
                name = "BP" + util::toString(counter);
                auto p2 = particleSystem->addParticle(name, catalog_->staticBP());
                if (plane == Plane::XY) {
                    position_t r{coord1, coord2, box->lengthZ()};
                    p2->position(r);
                } else if (plane == Plane::YZ) {
                    position_t r{box->lengthX(), coord1, coord2};
                    p2->position(r);
                } else {
                    position_t r{coord2, box->lengthY(), coord1};
                    p2->position(r);
                }
                particles.emplace_back(p2);
                coord2 += dis2();
            } while (coord2 <= limit2());
            coord2 = excludeCorner ? dis2() : 0.0;
            coord1 += dis1();
        } while (coord1 <= limit1());
        logger.info("Added number of boundary particles: " + util::toString(counter));

        particleSystem->addParticleGroup(particles, bonds);
    }

    p_system_ptr_t
    ParticleSystemFactory::fromPDB(std::string &fileName,
                                   const temperature_t& temperature,
                                   bool excludeWater) {
        util::Logger logger("Creating atomic particle system from PDB entry.");
        auto source = std::make_shared<InputSource>(fileName);
        p_system_ptr_t particleSystem = factory::atomistic();
        cont_handler_ptr_t handler(new AtomicContentHandler(particleSystem,
                                                                    catalog_,
                                                                            excludeWater));
        PDBReader reader;
        reader.parse(handler, source);
        particleSystem->doWithAll<void>([temperature] (const std::vector<p_ptr_t>& all) {
            for (auto p : all) {
                util::assignVelocity(p, temperature);
            }
        });
        return particleSystem;
    }

    spec_catalog_ptr_t&
    ParticleSystemFactory::catalog() {
        return catalog_;
    }

}