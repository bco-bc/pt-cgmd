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
#include "simploce/util/util.hpp"
#include "simploce/units/units-mu.hpp"
#include "simploce/util/logger.hpp"
#include <stdexcept>
#include <cmath>
#include <utility>

namespace simploce {

    ParticleSystemFactory::ParticleSystemFactory(spec_catalog_ptr_t  catalog) :
        catalog_{std::move(catalog)}
    {
        if ( !catalog_ ) {
            util::Logger logger("simploce::ParticleSystemFactory::ParticleSystemFactory");
            util::logAndThrow(logger,
                              "simploce::ParticleSystemFactory: Missing particle specifications catalog.");
        }
    }

    at_sys_ptr_t
    ParticleSystemFactory::diatomic(const distance_t& distance,
                                    const spec_ptr_t& spec,
                                    const temperature_t& temperature) {
        util::Logger logger{"simploce::ParticleSystemFactory::diatomic"};

        auto atomistic = factory::atomistic();
        std::vector<p_ptr_t> atoms{};
        position_t r1{-0.5 * distance(), 0.0, 0.0};
        std::string name1 = spec->name() + "1";
        auto atom_1 = atomistic->addAtom(name1, spec);
        atom_1->position(r1);
        util::assignVelocity(atom_1, temperature);
        atoms.emplace_back(atom_1);
        position_t r2{0.5 * distance(), 0.0, 0.0};
        std::string name2 = spec->name() + "2";
        auto atom_2 = atomistic->addAtom(name2, spec);
        atom_2->position(r2);
        util::assignVelocity(atom_2, temperature);
        atoms.emplace_back(atom_2);
        id_pair_t bond = std::make_pair(atom_1->id(), atom_2->id());
        std::vector<id_pair_t> bonds{bond};
        auto atomGroup = ParticleGroup::make(atoms, bonds);
        atomistic->addParticleGroup(atomGroup);
        box_ptr_t box(new Cube<real_t>{0.0});
        atomistic->box(box);
        return std::move(atomistic);
    }

    cg_sys_ptr_t
    ParticleSystemFactory::polarizableWater(const box_ptr_t& box,
                                            const density_t& densitySI,
                                            const temperature_t& temperature,
                                            std::size_t nLimit) {
        const distance_t DISTANCE_CW_DP = 0.2; // nm.

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
        ncgWaters = int(ncgWaters > nLimit ? nLimit : ncgWaters);
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
        const length_t spacing{0.53};  // Roughly the location of the first peak of g(r), in nm.
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
                    int l = int(util::random<real_t>() * 3.0);
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
        return std::move(coarseGrained);
    }

    at_sys_ptr_t
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
        const length_t spacing = 0.4; // The location of the first peak of g(r) for
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
                    real_t x = (i + 0.5) * spacing() + util::random<real_t>() * 0.1;
                    real_t y = (j + 0.5) * spacing() - util::random<real_t>() * 0.1;
                    real_t z = (k + 0.5) * spacing() + util::random<real_t>() * 0.1;
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
        return std::move(atomistic);
    }

    at_sys_ptr_t
    ParticleSystemFactory::argon(const box_ptr_t& box,
                                 const density_t& densitySI,
                                 const temperature_t& temperature)
    {
        util::Logger logger{"ParticleSystemFactory:::argon"};

        logger.info("Creating atomistic particle model for liquid argon.");
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
        logger.debug("Spacing between LJ beads: " + util::toString(spacing));

        // Box dimensions.
        length_t Lx = box->lengthX();
        length_t Ly = box->lengthY();
        length_t Lz = box->lengthZ();
        volume_t volume = box->volume();
        logger.info("Requested box size (nm): " + util::toString(box->size()));
        logger.debug("Box volume (nm^3): " + util::toString(volume));

        int nAtoms = int(numberDensity() * volume());
        logger.info("Requested number of argon atoms: " + util::toString(nAtoms));

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
        logger.info("Density (u/nm3): " + util::toString(spec->mass() * atomistic->numberOfAtoms() / volume()));
        return std::move(atomistic);
    }

    spec_catalog_ptr_t&
    ParticleSystemFactory::catalog() {
        return catalog_;
    }

}