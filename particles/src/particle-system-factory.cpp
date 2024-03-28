/*
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 12:50 PM
 */

#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/p-util.hpp"
#include "simploce/particle/atomic-content-handler.hpp"
#include "simploce/conf/p-conf.hpp"
#include "simploce/chem/input-source.hpp"
#include "simploce/chem/pdb-reader.hpp"
#include "simploce/util/util.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include <stdexcept>
#include <cmath>
#include <utility>
#include <random>
#include <vector>
#include <algorithm>

namespace simploce {
    namespace particle_factory {

        /**
         * Attempts to neutralize given particle system by adding positive or negative ions.
         * @param particleSystem Particle system.
         * @param catalog Particle specification catalog.
         * @return True if successful, otherwise false.
         */
        static
        bool neutralize(const p_system_ptr_t& particleSystem,
                        const spec_catalog_ptr_t& catalog) {
            util::Logger logger("simploce::particle_factory::neutralize()");
            logger.trace("Entering.");

            auto box = particleSystem->box();
            bool success{true};
            auto totalCharge = particleSystem->charge();
            if (totalCharge() > 0.0) {
                int nAdd = int(totalCharge());
                auto Cl = catalog->lookup("Cl-");
                logger.debug(Cl->name() + ": Adding this specification.");
                for (int n = 0; n != nAdd; ++n) {
                    auto p = particleSystem->addParticle("Cl-", Cl);
                    auto x = util::randomUniform() * box->lengthX();
                    auto y = util::randomUniform() * box->lengthY();
                    auto z = util::randomUniform() * box->lengthZ();
                    p->position(position_t{x, y, z});
                }
                logger.info(std::to_string(nAdd) + ": Number of extra Cl- ions added.");
            } else if (totalCharge() < 0.0) {
                int nAdd = int(std::fabs(totalCharge()));
                auto Na = catalog->lookup("Na+");
                logger.debug(Na->name() + ": Adding this specification.");
                for (int n = 0; n != nAdd; ++n) {
                    auto p = particleSystem->addParticle("Na+", Na);
                    auto x = util::randomUniform() * box->lengthX();
                    auto y = util::randomUniform() * box->lengthY();
                    auto z = util::randomUniform() * box->lengthZ();
                    p->position(position_t{x, y, z});
                }
                logger.info(std::to_string(nAdd) + ": Number of extra Na+ ions added.");
            }

            totalCharge = particleSystem->charge();
            if (std::fabs(totalCharge()) > 0.0) {
                logger.warn(util::to_string(totalCharge) + ": Still non-zero total charge.");
                success = false;
            }
            logger.trace("Leaving.");
            return success;
        }

    }

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

        // Atomistic model
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
    ParticleSystemFactory::cgmdPolarizableWater(const box_ptr_t& box,
                                                std::size_t nLimit,
                                                const density_t& densitySI,
                                                const temperature_t& temperature) {
        const dist_t DISTANCE_CW_DP = 0.2; // nm.

        util::Logger logger{"ParticleSystemFactory::cgmdPolarizableWater"};

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
        logger.info("Generating " + std::to_string(ncgWaters) + " CG water.");
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
        logger.debug("Can create a maximum of " + std::to_string(nx * ny * nz) + " CW beads.");

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
        logger.info("Created " + std::to_string(coarseGrained->numberOfBeads()) + " beads.");
        logger.info("Created " + std::to_string(coarseGrained->numberOfParticleGroups()) + " water groups.");
        if ( coarseGrained->numberOfParticleGroups() < ncgWaters ) {
            logger.warn("The number of created water groups is less then the number of requested water groups.");
        }
        coarseGrained->box(box);

        util::removeOverallLinearMomentum(coarseGrained);

        return std::move(coarseGrained);
    }

    p_system_ptr_t
    ParticleSystemFactory::mesoscalePolarizableWater(const box_ptr_t &box,
                                                     int numberOfWaterFluidElements,
                                                     const temperature_t T) {
        util::Logger logger{"simploce::ParticleSystemFactory::mesoscalePolarizableWater"};

        logger.info(std::to_string(box->lengthX()) + ", " +
                     std::to_string(box->lengthY()) + ", " +
                     std::to_string(box->lengthZ()) + ": Box dimensions");
        logger.info(std::to_string(numberOfWaterFluidElements) + ": Requested number of water particles.");
        logger.info(std::to_string(T()) + ": Requested temperature.");
        auto volume = box->volume();
        density_t rho = real_t(2.0 * numberOfWaterFluidElements) / volume;
        logger.info(std::to_string(rho()) + ": Requested total bead number density.");

        auto particleSystem = factory::coarseGrained();
        auto CW = catalog_->lookup("CW");
        auto DP = catalog_->lookup("DP");

        // Position of DP drawn from normal distribution
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<real_t> d{0.0, 0.2};
        for (auto i = 0; i != numberOfWaterFluidElements; ++i) {
            // Single water group.
            std::vector<p_ptr_t> beads{};
            std::vector<id_pair_t> bonds{};

            auto x = util::randomUniform(0.0, box->lengthX());
            auto y = util::randomUniform(0.0, box->lengthY());
            auto z = util::randomUniform(0.0, box->lengthZ());
            auto name = "CW" + std::to_string(i);
            auto cwBead = particleSystem->addBead(name, CW);
            cwBead->position(position_t{x, y, z});
            util::assignVelocity(cwBead, T, true);
            beads.emplace_back(cwBead);

            name = "DP" + std::to_string(i);
            auto dpBead = particleSystem->addBead(name, DP);
            x += d(gen);
            y += d(gen);
            z += d(gen);
            dpBead->position(position_t{x, y, z});
            util::assignVelocity(dpBead, T, true);
            beads.emplace_back(dpBead);

            id_pair_t pair = std::make_pair(cwBead->id(), dpBead->id());
            bonds.emplace_back(pair);
            particleSystem->addParticleGroup(beads, bonds);
        }
        particleSystem->box(box);

        util::removeOverallLinearMomentum(particleSystem);

        return std::move(particleSystem);
    }

    p_system_ptr_t
    ParticleSystemFactory::simpleElectrolyte(const box_ptr_t& box,
                                             const molarity_t& molarity,
                                             const temperature_t& temperature,
                                             bool randomPositions)
    {
        util::Logger logger{"ParticleSystemFactory::simpleElectrolyte"};
        logger.trace("Entering");

        logger.info("Creating electrolyte particle model for a simple electrolyte solution:");

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
        if (!randomPositions) {
            logger.debug("Distance spacing between ion particles: " + boost::lexical_cast<std::string>(spacing));
            logger.debug("Number of coordinates in x-direction: " + boost::lexical_cast<std::string>(nx));
            logger.debug("Number of coordinates in y-direction: " + boost::lexical_cast<std::string>(ny));
            logger.debug("Number of coordinates in z-direction: " + boost::lexical_cast<std::string>(nz));
        }

        // Start from an empty particle model.
       auto electrolyte = factory::atomistic();

        // Create ions.
        std::size_t counter = 0;
        int i = 0, j = 0, k = 0;

        while ( i < nx && counter < nIons ) {
            while ( j < ny && counter < nIons ) {
                while ( k < nz && counter < nIons ) {

                    // Position.
                    real_t x, y, z;
                    if (randomPositions) {
                        x = util::randomUniform() * Lx();
                        y = util::randomUniform() * Ly();
                        z = util::randomUniform() * Lz();
                    } else {
                        x = (i + 0.5) * spacing() + util::random() * 0.1;
                        y = (j + 0.5) * spacing() - util::random() * 0.1;
                        z = (k + 0.5) * spacing() + util::random() * 0.1;
                    }
                    position_t r{x,y,z};
                    if ( counter % 2 == 0 ) {
                        // Na+
                        std::string name = NA->name() + std::to_string(counter + 1 );
                        auto ion = electrolyte->addAtom(name, NA);
                        ion->position(r);
                        util::assignVelocity(ion, temperature);
                    } else {
                        // Cl-
                        std::string name = CL->name() + std::to_string(counter + 1);
                        auto ion = electrolyte->addAtom(name, CL);
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
        electrolyte->box(box);
        particle_factory::neutralize(electrolyte, catalog_);

        logger.info("Number of ion particles generated: " + boost::lexical_cast<std::string>(counter));
        if (counter < nIons) {
            logger.warn("Number of generated ions is smaller than requested.");
        }
        logger.info("TOTAL ion number density (1/nm^3): " + boost::lexical_cast<std::string>(real_t(counter)/volume()));

        util::removeOverallLinearMomentum(electrolyte);

        return std::move(electrolyte);
    }

    p_system_ptr_t
    ParticleSystemFactory::largeObjectInElectrolyte(const box_ptr_t& box,
                                                    const molarity_t& molarity,
                                                    const temperature_t& temperature,
                                                    const std::string &specName) {
        util::Logger logger{"simploce::ParticleSystemFactory::largeObjectInElectrolyte"};
        logger.trace("Entering");

        logger.info("Creating atomistic particle model for a simple electrolyte solution "
                    "with a large object at the center.");

        logger.info("Molarity (mol/l): " + boost::lexical_cast<std::string>(molarity()));
        logger.info("Temperature (K) : " + boost::lexical_cast<std::string>(temperature()));
        logger.info(specName + ": Specification name of large object.");

        // Create electrolyte solution.
        auto particleSystem = this->simpleElectrolyte(box, molarity, temperature);

        // Place large object at the center of the box.
        auto center = particleSystem->box()->center();
        auto spec = catalog_->lookup(specName);
        auto particle = particleSystem->addParticle(specName, spec);
        particle->position(center);

        // Neutralize.
        particle_factory::neutralize(particleSystem, catalog_);
        logger.info(std::to_string(particleSystem->charge()()) + ": Total charge.");

        logger.info(std::to_string(particleSystem->numberOfParticles()) + ": Total number of particles.");

        // Done
        logger.trace("Leaving");
        return particleSystem;
    }

    p_system_ptr_t
    ParticleSystemFactory::argon(const box_ptr_t& box,
                                 std::size_t nLimit,
                                 const density_t& densitySI,
                                 const temperature_t& temperature)
    {
        util::Logger logger{"simploce::ParticleSystemFactory::argon"};

        logger.info("Creating atomistic particle model for argon.");
        logger.info("Requested temperature: " + std::to_string(temperature()));

        // Conversion of kg/m^3 to u/nm^3.
        density_t density = densitySI / (units::si<real_t>::MU * 1.0e+27);

        // Argon specification.
        auto spec = catalog_->lookup("Ar");
        logger.info("Argon specification:");
        util::log(logger, *spec);

        number_density_t numberDensity = density / spec->mass();
        logger.info("Requested density (kg/m^3): " + std::to_string(densitySI()));
        logger.info("Requested density (u/nm^3): " + std::to_string(density()));
        logger.info("Requested number density (1/nm^3):" + std::to_string(numberDensity()));

        // Spacing between argon atoms.
        const length_t spacing{0.3};
        logger.debug("Spacing between LJ atoms: " + std::to_string(spacing()));

        // Box dimensions.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        logger.info("Requested box size (nm): " + std::to_string(box->size()));
        logger.debug("Box volume (nm^3): " + std::to_string(volume()));

        int nAtoms = int(numberDensity() * volume());
        logger.info("Requested number of argon atoms: " + std::to_string(nAtoms));
        if (nAtoms > int(nLimit)) {
            logger.warn("Number of requested atoms is higher that maximum number of atoms allowed.");
            nAtoms = int(nLimit);
        }
        logger.info("Generating " + std::to_string(nAtoms) + " argon atoms.");

        int nx = util::nint(Lx / spacing);
        int ny = util::nint(Ly / spacing);
        int nz = util::nint(Lz / spacing);
        logger.debug("Number of coordinates in x-direction: " + std::to_string(nx));
        logger.debug("Number of coordinates in y-direction: " + std::to_string(ny));
        logger.debug("Number of coordinates in z-direction: " + std::to_string(nz));

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
                    std::string name = spec->name() + std::to_string(counter + 1);
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
        logger.info("Number of Argon atoms generated: " + std::to_string(atomistic->numberOfAtoms()));
        numberDensity = real_t(atomistic->numberOfAtoms()) / volume();
        logger.info("Number density (1/nm^3): " + std::to_string(numberDensity()));
        density = spec->mass()() * numberDensity();
        logger.info("Density (u/nm3): " + std::to_string(density()));
        density_t densSI = density() * units::si<real_t>::MU * 1.0e+27;
        logger.info("Density (kg/m^3): " + std::to_string(densSI()));

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
        logger.info(util::to_string(*box) + ": Requested box dimensions.");
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
        logger.info(std::to_string(placeRandom) + ": Place beads randomly in box?");

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
                    auto name = spec->name() + std::to_string(counter);
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
                    bonds.emplace_back(id1, id2);
                }
                polymerSystem->addParticleGroup(particles, bonds);
            }

            // Add water beads, placed randomly in the box.
            for (auto i = 0; i != numberOfWaters; ++i) {
                auto x = dis(gen) * boxX;
                auto y = dis(gen) * boxY;
                auto z = dis(gen) * boxZ;
                position_t r{x, y, z};
                auto name = waterBeadSpecName + std::to_string(counter);
                auto particle = polymerSystem->addBead(name, waterSpec);
                particle->position(r);
                util::assignVelocity(particle, temperature, true);
                counter += 1;
                if ( i > 0 && i % 10000 == 0) {
                    logger.debug(std::to_string(i) + ": Number of water beads created.");
                }
            }
            logger.debug(std::to_string(numberOfWaters) + ": Number of water beads created.");
            logger.info(std::to_string(counter) + ": Total number of beads placed randomly in box.");
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
                    auto name = spec->name() + std::to_string(counter);
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
                    bonds.emplace_back(id1, id2);
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
                auto name = waterBeadSpecName + std::to_string(counter);
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

    p_system_ptr_t
    ParticleSystemFactory::dropletPolymerSolution(const box_ptr_t& box,
                                                  int chainLength,
                                                  const std::string& monomericUnitBeadSpecName,
                                                  int numberOfPolymers,
                                                  const length_t& spacing,
                                                  int numberOfWaters,
                                                  const std::string& waterBeadSpecName,
                                                  int numberOfDropletBeads,
                                                  const std::string& dropletBeadSpecName,
                                                  const temperature_t& temperature) {
        util::Logger logger("simploce::ParticleSystemFactory::dropletPolymerSolution");
        logger.trace("Entering.");

        logger.info("Creating droplet polymer solution.");
        logger.info(util::to_string(*box) + ": Requested box dimensions.");
        logger.info(std::to_string(chainLength) + ": Requested polymer chain length (number of monomeric units).");
        logger.info(monomericUnitBeadSpecName + ": Monomeric unit particle specification name.");
        logger.info(std::to_string(numberOfPolymers) + ": Requested number of polymers.");
        logger.info(std::to_string(spacing()) + ": Requested spacing between monomeric units.");
        logger.info(std::to_string(numberOfWaters) + ": Requested number of water beads.");
        logger.info(waterBeadSpecName + ": Water specification name.");
        logger.info(std::to_string(temperature()) + ": Requested temperature.");
        logger.info(std::to_string(numberOfPolymers * chainLength) + ": Number of polymer beads required.");
        logger.info(std::to_string(numberOfDropletBeads) + ": Requested number of droplet beads.");
        logger.info(std::to_string(numberOfWaters) + ": Requested number of water beads.");
        auto numberOfBeads = numberOfWaters + numberOfPolymers * chainLength + numberOfDropletBeads;
        logger.info(std::to_string(numberOfBeads) + ": Total number of beads required.");

        auto particleSystem = this->polymerSolution(box,
                                                                 chainLength,
                                                                 monomericUnitBeadSpecName,
                                                                 numberOfPolymers,
                                                                 spacing,
                                                                 numberOfWaters,
                                                                 waterBeadSpecName,
                                                                 temperature,
                                                                 true);

        auto boxX = box->lengthX();
        auto boxY = box->lengthY();
        auto boxZ = box->lengthZ();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<real_t> dis(0.0, 1.0);
        auto dropletBeadSpec = catalog_->lookup(dropletBeadSpecName);
        auto counter = particleSystem->numberOfParticles();

        // Add droplet particles, placed randomly in the box.
        for (auto i = 0; i != numberOfDropletBeads; ++i) {
            auto x = dis(gen) * boxX;
            auto y = dis(gen) * boxY;
            auto z = dis(gen) * boxZ;
            position_t r{x, y, z};
            auto name = dropletBeadSpecName + std::to_string(counter);
            auto particle = particleSystem->addParticle(name, dropletBeadSpec);
            particle->position(r);
            util::assignVelocity(particle, temperature, true);
            counter += 1;
            if ( i > 0 && i % 10000 == 0) {
                logger.debug(std::to_string(i) + ": Number of droplet beads created.");
            }
        }
        logger.debug(std::to_string(numberOfDropletBeads) + ": Number of droplet beads created.");

        logger.trace("Leaving.");
        return std::move(particleSystem);
    }

    p_system_ptr_t
    ParticleSystemFactory::identicalParticles(const simploce::box_ptr_t &box,
                                              const std::string &specName,
                                              const simploce::number_density_t &rho,
                                              const temperature_t& temperature,
                                              bool mesoscale) {
        util::Logger logger("simploce::ParticleSystemFactory::identicalParticles()");
        logger.trace("Entering");

        logger.info("Creating physical system of identical particles.");
        logger.info(util::to_string(*box) + ": Requested box dimensions.");
        logger.info(specName + ": Requested particle specification name.");
        logger.info(std::to_string(rho()) + ": Requested number density.");
        logger.info(std::to_string(temperature()) + ": Requested temperature.");
        auto numberOfParticles = int(rho() * box->volume());
        logger.info(std::to_string(numberOfParticles) + ": Requested number of particles.");
        logger.info(std::to_string(mesoscale) + ": Mesoscale?");

        auto spec = catalog_->lookup(specName);

        auto particleSystem = factory::coarseGrained();
        for (auto counter = 0; counter != numberOfParticles; ++counter) {
            auto x = util::randomUniform(0.0, box->lengthX());
            auto y = util::randomUniform(0.0, box->lengthY());
            auto z = util::randomUniform(0.0, box->lengthZ());
            position_t r(x, y, z);
            std::string name = specName + std::to_string(counter);
            auto particle = particleSystem->addBead(name, spec);
            particle->position(r);
            util::assignVelocity(particle, temperature, mesoscale);
            auto n = particleSystem->numberOfParticles();
            if (n % 100000 == 0) {
                logger.info(std::to_string(n) + ": Number of particles created.");
            }
        }
        logger.info(std::to_string(particleSystem->numberOfParticles()) + ": Number of particles created.");
        particleSystem->box(box);

        logger.trace("Leaving.");
        return particleSystem;
    }

    void
    ParticleSystemFactory::addBoundaryParticles(const p_system_ptr_t& particleSystem,
                                                dist_t spacing,
                                                const spec_ptr_t& specBP,
                                                srf_charge_density_t sigma,
                                                Plane plane,
                                                bool bothSides,
                                                bool excludeCorner,
                                                temperature_t temperature,
                                                bool mesoscale,
                                                bool rough,
                                                length_t boundaryWidth) {
        util::Logger logger{"simploce::ParticleSystemFactory::addBoundaryParticles"};
        logger.info("Spacing: " + std::to_string(spacing()));
        logger.info(util::to_string(sigma) + ": Requested surface charge density.");
        logger.info("Plane: " + plane.toString());
        logger.info(std::to_string(bothSides) + ": Both sides of box?");
        logger.info(std::to_string(excludeCorner) + ": Exclude corner?");
        logger.info(util::to_string(temperature) + ": Temperature.");
        logger.info(std::to_string(mesoscale) + ": Mesoscale system?");
        logger.info(std::to_string(rough) + ": Rough surface.");
        logger.info(util::to_string(boundaryWidth) + ": Boundary width.");
        logger.info(specBP->name() + ": Boundary particle specification name.");

        auto box = particleSystem->box();
        length_t limit1, limit2, dis1, dis2;
        length_t displacement = std::sqrt(boundaryWidth() * boundaryWidth() / 3.0);
        logger.debug(std::to_string(displacement()) + ": Range for displacement of boundary particles.");

        int n1, n2; // Number of points in each direction.
        area_t area;
        if (plane == Plane::XY) {
            n1 = int(box->lengthX() / spacing());
            dis1 = box->lengthX() / n1;
            n2 = int(box->lengthY() / spacing());
            dis2 = box->lengthY() / n2;
            limit1 = box->lengthX();
            limit2 = box->lengthY();
            area = box->lengthX() * box->lengthY();
        } else if (plane == Plane::YZ) {
            n1 = int(box->lengthY() / spacing());
            dis1 = box->lengthY() / n1;
            n2 = int(box->lengthZ() / spacing());
            dis2 = box->lengthZ() / n2;
            limit1 = box->lengthY();
            limit2 = box->lengthZ();
            area = box->lengthY() * box->lengthZ();
        } else {
            n1 = int(box->lengthZ() / spacing());
            dis1 = box->lengthZ() / n1;
            n2 = int(box->lengthX() / spacing());
            dis2 = box->lengthX() / n2;
            limit1 = box->lengthZ();
            limit2 = box->lengthX();
            area = box->lengthZ() * box->lengthX();
        }
        logger.debug(std::to_string(n1) + ": Value of n1.");
        logger.debug(std::to_string(n2) + ": Value of n2.");
        logger.debug(util::to_string(area) + ": Surface area.");

        real_t coord1 = excludeCorner ? spacing() : 0.0;
        limit1 = excludeCorner ? limit1 - dis1 : limit1;
        real_t coord2 = excludeCorner ? spacing() : 0.0;
        limit2 = excludeCorner ? limit2 - dis2 : limit2;
        logger.debug(std::to_string(coord1) + ": Lower limit of coord1.");
        logger.debug(std::to_string(limit1()) + ": Upper limit for coord1.");
        logger.debug(std::to_string(coord2) + ": Lower limit of coord2.");
        logger.debug(std::to_string(limit2()) + ": Upper limit for coord2.");

        std::vector<p_ptr_t> particles;
        std::string nameBP = specBP->name();
        int counter = 0;
        do {
            do {
                counter += 1;

                // One side of box.
                std::string name = nameBP + std::to_string(counter);
                auto p1 = particleSystem->addParticle(name, specBP);
                if (plane == Plane::XY) {
                    position_t r1{coord1, coord2, 0.0};
                    p1->position(r1);
                } else if (plane == Plane::YZ) {
                    position_t r1{0.0, coord1, coord2};
                    p1->position(r1);
                } else {
                    position_t r1{coord2, 0.0, coord1};
                    p1->position(r1);
                }

                // Displace boundary particle?
                if (rough) {
                    auto x = util::randomUniform(-displacement(), displacement());
                    auto y = util::randomUniform(-displacement(), displacement());
                    auto z = util::randomUniform(-displacement(), displacement());
                    auto r = p1->position() + position_t{x, y, z};
                    p1->position(r);
                }

                // Velocity.
                util::assignVelocity(p1, temperature, mesoscale);
                particles.emplace_back(p1);

                // Update counter.
                counter += 1;

                // Other side of box.
                if (bothSides) {
                    name = "BP" + std::to_string(counter);
                    auto p2 = particleSystem->addParticle(name, specBP);
                    if (plane == Plane::XY) {
                        position_t r2{coord1, coord2, box->lengthZ()};
                        p2->position(r2);
                    } else if (plane == Plane::YZ) {
                        position_t r2{box->lengthX(), coord1, coord2};
                        p2->position(r2);
                    } else {
                        position_t r2{coord2, box->lengthY(), coord1};
                        p2->position(r2);
                    }

                    // Displace boundary particle?
                    if (rough) {
                        auto x = util::randomUniform(-displacement(), displacement());
                        auto y = util::randomUniform(-displacement(), displacement());
                        auto z = util::randomUniform(-displacement(), displacement());
                        auto r = p2->position() + position_t{x, y, z};
                        p2->position(r);
                    }

                    // Velocity.
                    util::assignVelocity(p2, temperature, mesoscale);
                    particles.emplace_back(p2);

                     // Update counter.
                    counter += 1;
                }
                coord2 += dis2();
            } while (coord2 <= limit2());
            coord2 = excludeCorner ? spacing() : 0.0;
            coord1 += dis1();
        } while (coord1 <= limit1());

        auto boundaryParticles = particleSystem->ofSpec(specBP);
        auto numberOfBP = boundaryParticles.size();
        for (auto &p: boundaryParticles) {
            particleSystem->resetToFree(p);
        }
        logger.info(std::to_string(numberOfBP) + ": Number of boundary particles generated.");
        charge_t chargeBP = area() * sigma() / real_t(numberOfBP);
        logger.debug(util::to_string(chargeBP) + ": Charge per boundary particle.");

        // Replace boundary particle specification to reflected requested surface density.
        auto spec = ParticleSpec::create(specBP->name(),
                                         chargeBP,
                                         specBP->mass(),
                                         specBP->radius(),
                                         specBP->isFree(),
                                         specBP->description());
        logger.debug(util::to_string(*spec) + ": Specification newly created surface boundary particle.");
        logger.trace("Resetting specification of boundary particle...");
        for (auto& p: boundaryParticles) {
            particleSystem->resetSpec(p, spec);
        }
        logger.trace("Done");

        // Neutralize.
        auto success = particle_factory::neutralize(particleSystem, catalog_);
        if (!success) {
            ParticleSystemFactory::adjustBoundaryParticleToNeutralize(particleSystem, spec);
        }

        // Adjust momenta.
        util::removeOverallLinearMomentum(particleSystem);

        // Done.
        logger.trace("Leaving");
    }

    srf_charge_density_t
    ParticleSystemFactory::addChargedSurface(const p_system_ptr_t& particleSystem,
                                             srf_charge_density_t sigma,
                                             temperature_t temperature) {
        util::Logger logger("simploce::ParticleSystemFactory::addChargedSurface()");
        logger.trace("Entering");

        logger.debug(util::to_string(sigma) + ": Requested surface density.");
        logger.debug(util::to_string(temperature) + ": Requested temperature.");

        auto box = particleSystem->box();
        auto area = box->lengthX() * box->lengthY();
        auto charge = util::nint(area * sigma());
        logger.debug(util::to_string(charge) + ": Charge 'fake' surface particle.");
        auto spec = ParticleSpec::create("SBP", charge, 1.0, 1.0, true, "SBP");
        auto sbp = particleSystem->addParticle("SBP", spec);
        auto success = particle_factory::neutralize(particleSystem, catalog_);
        if (!success) {
            util::logAndThrow(logger, "Cannot neutralize given particle system.");
        }
        particleSystem->remove(sbp);  // No need to keep.
        srf_charge_density_t adjustedSigma = real_t(charge) / area;
        logger.debug(util::to_string(adjustedSigma) + ": Adjusted surface charge density.");

        logger.trace("Leaving");
        return adjustedSigma;
    }

    void
    ParticleSystemFactory::makeChannel(const simploce::p_system_ptr_t &particleSystem,
                                       const simploce::length_t &wallWidth,
                                       bool mesoscale,
                                       const number_density_t& rho,
                                       bool adjustToNumberDensity,
                                       const temperature_t& temperature) {
        util::Logger logger("simploce::ParticleSystemFactory::makeChannel()");
        logger.trace("Entering.");
        logger.info(std::to_string(wallWidth()) + ": Requested wall width.");
        auto box = particleSystem->box();
        logger.info(util::to_string(*box) + ": Box dimensions.");
        logger.info(std::to_string(mesoscale) +  ": Mesoscale?");

        // Specification for a boundary particle.
        auto SBP = catalog_->staticBP();

        particleSystem->doWithAllFreeGroups<void>([box, wallWidth, particleSystem, SBP, mesoscale, temperature] (
                std::vector<p_ptr_t>& all,
                std::vector<p_ptr_t>& free,
                std::vector<pg_ptr_t>& groups) {
            util::Logger logger("simploce::ParticleSystemFactory::makeChannel()");
            logger.debug(std::to_string(all.size()) + ": Total number of particles.");
            // Particles in particle groups.
            logger.debug(std::to_string(groups.size()) + ": Number of particle groups.");
            std::size_t counter = 0;
            std::vector<p_ptr_t> makeFree{};
            std::vector<size_t> removeGroup{};
            for (auto& g : groups) {
                bool groupIsWall{false};
                std::vector<id_t> ids{};

                // Check individual particles.
                auto particles = g->particles();
                for (auto& p: particles) {
                    auto r = p->position();
                    auto x = std::fabs(r.x());
                    auto y = std::fabs(r.y());

                    // Should the given particle be converted to a boundary particle?
                    auto inWall = (box->lengthX() - x < wallWidth() ||
                                  x < wallWidth() ||
                                  box->lengthY() - y < wallWidth() ||
                                  y < wallWidth());
                    if (inWall) {
                        ids.emplace_back(p->id());
                        groupIsWall = true;  // Mark the group as being part of the wall.
                    }
                }

                // Only if the group is marked as being part of the boundary, convert all its particles to
                // boundary particles.
                if (groupIsWall) {

                    // The group ceases to exist as a group, because the group's particles are now individual
                    // boundary particles. Mark the group for removal.
                    removeGroup.emplace_back(g->id());

                    // Convert.
                    for (auto& p: particles) {
                        std::string name = "BP" + std::to_string(counter);
                        particleSystem->resetSpec(p, SBP);
                        particleSystem->resetName(p, name);
                        util::assignVelocity(p, temperature, mesoscale);
                        makeFree.emplace_back(p);
                        counter += 1;
                    }
                }
            }
            logger.info(std::to_string(counter) + ": Number of boundary particles from particle groups.");

            // Any free particles that should be converted to boundary particles.
            counter = 0;
            for (auto& p : free) {
                auto r = p->position();
                auto x = std::fabs(r.x());
                auto y = std::fabs(r.y());
                auto inWall = box->lengthX() - x < wallWidth() ||
                              x < wallWidth() ||
                              box->lengthY() - y < wallWidth() ||
                              y < wallWidth();
                if (inWall) {
                    std::string name = "BP" + std::to_string(counter);
                    particleSystem->resetSpec(p, SBP);
                    particleSystem->resetName(p, name);
                    util::assignVelocity(p, temperature, mesoscale);
                    makeFree.emplace_back(p);
                    counter += 1;
                }
            }
            logger.info(std::to_string(counter) + ": Number of boundary particles from free particles.");

            // Remove groups whose particles are now boundary particles. The latter become free particles.
            std::vector<pg_ptr_t> groupsToRemove{};
            for (auto& g : groups) {
                auto id = g->id();
                if (std::find(removeGroup.begin(), removeGroup.end(), id) != removeGroup.end()) {
                   groupsToRemove.emplace_back(g);
                }
            }
            for (auto& g: groupsToRemove) {
                particleSystem->removeFromFree(g);
                particleSystem->removeGroup(g);
            }
            for (auto& p: makeFree) {
                particleSystem->resetToFree(p);
            }
        });

        logger.debug(std::to_string(particleSystem->numberOfFreeParticles()) + ": Current number of free particles.");
        logger.debug(std::to_string(particleSystem->numberOfSpecifications(SBP)) + ": Current number of boundary particles.");
        if (adjustToNumberDensity)
            this->adjustNumberDensityByRemovingParticleGroups(particleSystem, rho);

        // Balance the total momentum.
        util::removeOverallLinearMomentum(particleSystem);

        logger.info(std::to_string(particleSystem->numberOfFreeParticles()) + ": Number of free particles.");
        logger.info(std::to_string(particleSystem->numberOfSpecifications(SBP)) + ": Number of boundary particles.");

        // Done.
        logger.trace("Leaving.");
    }

    p_system_ptr_t
    ParticleSystemFactory::fromPDB(std::string &fileName,
                                   const temperature_t& temperature,
                                   bool excludeWater) {
        util::Logger logger("Creating atomic particle system from PDB entry.");
        auto source = std::make_shared<InputSource>(fileName);
        p_system_ptr_t particleSystem = factory::atomistic();
        cont_handler_ptr_t handler(new AtomicContentHandler(particleSystem, catalog_, excludeWater));
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

    void
    ParticleSystemFactory::adjustNumberDensityByRemovingParticleGroups(const p_system_ptr_t& particleSystem,
                                                                       const number_density_t& rho) {
        util::Logger logger{"simploce::ParticleSystemFactory::adjustNumberDensityByRemovingParticleGroups()"};
        logger.trace("Entering");

        logger.debug(std::to_string(real_t(particleSystem->numberOfParticles()) / particleSystem->box()->volume()) +
                     ": Total number density.");
        logger.debug(std::to_string(particleSystem->box()->volume()) + ": Box volume.");
        auto specSBP = catalog_->staticBP();
        auto numberOfSBP = particleSystem->numberOfSpecifications(specSBP);
        logger.info(std::to_string(particleSystem->numberOfParticles() - numberOfSBP) +
                    ": Number of particles other than boundary particles.");
        auto diameter = 2.0 * specSBP->radius();
        auto box = *particleSystem->box();
        real_t volume = (box[0] - diameter()) * (box[1] - diameter()) * box[2];
        logger.debug(std::to_string(volume) + ": Effective box volume.");

        number_density_t current = real_t((particleSystem->numberOfParticles() - numberOfSBP)) / volume;
        logger.debug(std::to_string(current()) + ": Initial number density for other than boundary particles.");
        auto numberOfGroups = particleSystem->numberOfParticleGroups();
        while (current > rho && numberOfGroups > 0) {
            auto group = particleSystem->doWithAllFreeGroups<pg_ptr_t>([](
                std::vector<p_ptr_t> &all,
                std::vector<p_ptr_t> &free,
                std::vector<pg_ptr_t> &groups
            ) {
                return groups[0];
            });
            auto particles = group->particles();
            particleSystem->removeGroup(group);
            for (auto& p: particles) {
                particleSystem->remove(p);
            }
            numberOfGroups = particleSystem->numberOfParticleGroups();
            current = real_t((particleSystem->numberOfParticles() - numberOfSBP)) / volume;
            logger.debug(std::to_string(particleSystem->numberOfParticles()) + ": Total number of particles.");
            logger.debug(std::to_string(current()) + ": Current number density for other than boundary particles.");
        }
        logger.info(std::to_string(current()) + ": Effective final number density for other than boundary particles.");
        logger.info(std::to_string(particleSystem->numberOfParticles() - numberOfSBP) +
                    ": Final number of particles other than boundary particles.");

        logger.trace("Leaving");
    }

    /**
     * Adjust properties of the boundary particles to neutralize particle system.
     * @param particleSystem Particle system.
     * @param specBP Boundary particle specification
     */
    void
    ParticleSystemFactory::adjustBoundaryParticleToNeutralize(const p_system_ptr_t& particleSystem,
                                                              const spec_ptr_t& specBP) {
        util::Logger logger("simploce::particle_factory::adjustBoundaryParticleToNeutralize()");
        logger.trace("Entering.");

        auto charge = particleSystem->charge();
        logger.debug(util::to_string(charge) + ": Total charge -before- adjustment.");
        if ( std::fabs(charge()) > 0.0) {
            auto specs = particleSystem->specsInUse();
            auto iter = specs.find(specBP->name());
            if (iter == specs.end()) {
                logger.warn(specBP->name() + ": No such particle specification in particle system.");
                logger.warn("Cannot adjust boundary particles. Particle system may not be neutral.");
                return;
            }
            auto n = particleSystem->numberOfSpecifications(specBP);
            logger.debug(std::to_string(n) + ": Number of particles of specification '" + specBP->name() + "'.");
            logger.debug(util::to_string(*specBP) + ": CURRENT surface boundary particle specification.");
            charge_t chargePerBP = std::fabs(charge() / real_t(n));
            auto adjustment = charge() > 0.0 ? -chargePerBP() : chargePerBP;
            auto newCharge = specBP->charge() + adjustment;
            auto spec = ParticleSpec::create(specBP->name(),
                                             newCharge,
                                             specBP->mass(),
                                             specBP->radius(),
                                             specBP->isFree(),
                                             specBP->description());
            logger.debug(util::to_string(adjustment) + ": Charge adjustment per surface boundary particle.");
            logger.debug(util::to_string(*spec) + ": ADJUSTED surface boundary particle specification.");
            auto particles = particleSystem->ofSpec(specBP);
            for (auto& p: particles) {
                particleSystem->resetSpec(p, spec);
            }
        }
        charge = particleSystem->charge();
        logger.debug(util::to_string(charge) + ": Total charge -after- adjustment.");
    }


}