/*
 * Converts particle system to PDB.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on April 20, 2022, 2.48 pm.
 */

#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-group.hpp"
#include "simploce/particle/bond.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/util.hpp"
#include "simploce/util/hybrid_36_c.hpp"
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace po = boost::program_options;
using namespace simploce;

static void
outputFreeParticle(std::ofstream& ostream,
                   const spec_catalog_ptr_t& catalog,
                   const p_ptr_t& particle,
                   int rsn,
                   const PeriodicBoundaryCondition& pbc,
                   bool inBox) {
    ostream << "ATOM  ";
    auto index = particle->index() + 1;
    auto stringIndex = std::to_string(index);
    if (index > 99999) {
        char result[5];
        const char* errmsg = hy36encode(5, index, result);
        if (errmsg) {
            throw std::domain_error(errmsg);
        }
        stringIndex = std::string(result);
    }
    ostream << std::setw(5) << stringIndex;
    //ostream << std::setw(5) << index;  // Atom serial number.

    ostream << conf::SPACE;      // Space between atom serial number and atom name.

    std::string name = particle->spec()->name();
    auto l = name.length() > 4 ? 4 : name.length();  // Atom name.
    for (std::size_t k = 0; k != l; ++k) {
        ostream << name[k];
    }
    for (std::size_t k = l; k < 4; ++k) {
        ostream << " ";
    }

    ostream << conf::SPACE; // Alternate location indicator.

    l = name.length() > 3 ? 3 : name.length();   // Residue name
    for (std::size_t k = 0; k != l; ++k) {
        ostream << name[k];
    }
    for (std::size_t k = l; k < 3; ++k) {
        ostream << " ";
    }

    ostream << conf::SPACE;

    ostream << conf::SPACE;                     // Chain identifier.
    int number = (rsn >= 10000 ? 9999 : rsn);
    ostream << std::setw(4) << number;  // Residue sequence number.
    ostream << conf::SPACE;             // Code for insertion of residues.

    ostream << "   ";                    // Space between code and coordinates.

    auto r_out = particle->position();
    position_t r;
    if (name == catalog->staticBP()->name()) {
        r = r_out;
    } else {
        if (inBox) {
            r = pbc.placeInside(r_out);
        } else {
            r = r_out;
        }
    }
    //r *= 10.0;  // To Angstrom.
    ostream.precision(3);
    ostream << std::setw(8) << r[0] << std::setw(8) << r[1] << std::setw(8) << r[2];

    // Occupancy, tempFactor
    ostream.precision(2);
    ostream << std::setw(6) << 0.0 << std::setw(6) << 0.0;

    for (auto k = 0; k != 10; ++k) {
        ostream << conf::SPACE;
    }

    ostream << conf::SPACE;
    ostream << particle->spec()->name()[0];
    ostream << "0.0";
    ostream << std::endl;
}


/**
 * Writes PDB.
 * @param particleSystem Particle system.
 * @param fnOutputPDB  Output file name.
 * @param catalog Particle specifications catalog.
 */
static
void toPDB(const p_system_ptr_t& particleSystem,
           const std::string& fnOutputPDB,
           const spec_catalog_ptr_t& catalog,
           bool inBox) {
    box_ptr_t box = particleSystem->box();
    PeriodicBoundaryCondition pbc(box);

    std::ofstream ostream;
    util::open_output_file(ostream, fnOutputPDB);
    ostream.setf(std::ios::fixed);
    ostream.precision(3);
    ostream << "HEADER    ";
    ostream << "PARTICLE SYSTEM PDB                    ";
    ostream << "06-APR-22 ";
    ostream << "  1PDB" << std::endl;
    particleSystem->doWithAllFreeGroups<void>([&ostream, &pbc, &catalog, inBox] (
            const std::vector<p_ptr_t>& all,
            std::vector<p_ptr_t>& free,
            std::vector<pg_ptr_t>& groups) {
        for (std::size_t i = 0; i != groups.size(); ++i) {
            auto& group = groups[i];
            std::string residueName = "GRP";
            auto& particles = group->particles();
            for (auto& p: particles) {
                ostream << "ATOM  ";
                int asn = int(p->index() + 1);
                ostream << std::setw(5) << asn;
                ostream << "  ";
                std::string name = p->spec()->name();
                int l = int(name.length() > 3 ? 3 : name.length());
                for (int k = 0; k < l; ++k) {
                    ostream << name[k];
                }
                for (int k = l; k < 3; ++k) {
                    ostream << " ";
                }
                ostream << " ";
                ostream << residueName;                // Residue name.
                ostream << " ";
                ostream << "A";                      // Chain identifier.
                auto rsn = i + 1;        // Residue sequence number.
                ostream << std::setw(4) << rsn;
                ostream << "    ";
                auto r_out = p->position();
                position_t r;
                if (name == catalog->staticBP()->name()) {
                    r = r_out;
                } else {
                    if (inBox) {
                        r = pbc.placeInside(r_out);
                    } else {
                        r = r_out;
                    }
                }

                //r *= 10.0;  // to Angstrom.
                ostream.precision(3);
                ostream << std::setw(8) << r[0] << std::setw(8) << r[1] << std::setw(8) << r[2];
                ostream.precision(2);
                ostream << conf::SPACE << std::setw(5) << 0.0 << conf::SPACE << std::setw(5) << 0.0;
                ostream << "           " << p->spec()->name()[0];
                ostream << std::endl;
            }
        }

        // Write free particles.
        int rsn = groups.size();
        for (auto& p : free) {
            std::string residueName = p->spec()->name().substr(0,3);
            rsn += 1;
            //int rsn = int (p->index() + 1);
            outputFreeParticle(ostream, catalog, p, rsn, pbc, inBox);
            // oneLineToPDB(ostream, "HETATM", residueName, rsn, p, catalog, pbc, inBox);
        }
    });
    ostream << "END" << std::endl;
    ostream.close();
}

int main(int argc, char *argv[]) {
    util::Logger logger{"simploce::s-pdb::main"};
    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnInputParticleSystem{"in.ps"};
    std::string fnOutputPDB{"out.pdb"};
    bool isCoarseGrained{false};
    bool inBox{true};

    po::options_description usage("Usage");
    usage.add_options() (
        "fn-particle-spec-catalog,s",
        po::value<std::string>(&fnParticleSpecCatalog),
        "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
    )(
        "fn-input-particle-system,i",
        po::value<std::string>(&fnInputParticleSystem),
        "Input file name particle system. Default is 'in.ps'."
    )(
        "coarse-grained,c",
        "Input is a coarse-grained description."
    )(
        "fn-output-pdb,o",
        po::value<std::string>(&fnOutputPDB),
        "Output file name PDB."
    )(
            "do-not-place-inside-box",
            "Do not place particles back into the central box."
    )(
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
        std::cout << "Convert particle system to PDB" << std::endl;
        std::cout << usage << "\n";
        return 0;
    }
    if (vm.count("fn-particle-spec-catalog")) {
        fnParticleSpecCatalog = vm["fn-particle-spec-catalog"].as<std::string>();
    }
    if (vm.count("fn-input-particle-system")) {
        fnInputParticleSystem = vm["fn-input-particle-system"].as<std::string>();
    }
    if (vm.count("coarse-grained") ) {
        isCoarseGrained = true;
    }
    if (vm.count("fn-output-pdb")) {
        fnOutputPDB = vm["fn-output-pdb"].as<std::string>();
    }
    if (vm.count("do-not-place-inside-box")) {
        inBox = false;
    }
    if (vm.count("verbose") ) {
        util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    }

    // Read particle specifications.
    auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);

    // Read particle system.
    auto particleSystem = factory::particleSystem(fnInputParticleSystem, catalog, isCoarseGrained);
    logger.info("Read particle system from '" + fnInputParticleSystem + "'.");
    logger.info("Number of particles: " + util::toString(particleSystem->numberOfParticles()));
    logger.info("Number of free particles: " + util::toString(particleSystem->numberOfFreeParticles()));
    logger.info("Number of particle groups: " + util::toString(particleSystem->numberOfParticleGroups()));

    // Write PDB.
    toPDB(particleSystem, fnOutputPDB, catalog, inBox);
    logger.info("Wrote output PDB to: " + fnOutputPDB);

    return EXIT_SUCCESS;
}