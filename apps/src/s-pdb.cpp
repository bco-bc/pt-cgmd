/*
 * Converts particle system to PDB.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on April 20, 2022, 2.48 pm.
 */

#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/pbc.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/util.hpp"
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace po = boost::program_options;
using namespace simploce;

void toPDB(const p_system_ptr_t& particleSystem,
           const std::string& fnOutputPDB,
           const spec_catalog_ptr_t& catalog) {
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
    particleSystem->doWithAllFreeGroups<void>([&ostream, &pbc, &catalog] (
            const std::vector<p_ptr_t>& all,
            std::vector<p_ptr_t>& free,
            std::vector<pg_ptr_t>& groups) {
        int counter = 1;
        for (auto p: all) {
            ostream << "HETATM";
            int asn = p->index() + 1;
            ostream << std::setw(5) << asn;
            ostream << " ";
            std::string name = p->spec()->name();
            int l = name.length() > 4 ? 4 : name.length();
            for (int k = 0; k != l; ++k) {
                ostream << name[k];
            }
            for (int k = l; k!= 4; ++k) {
                ostream << " ";
            }
            ostream << " ";
            ostream << "RES";                    // Residue name.
            ostream << " ";
            ostream << " ";                      // Chain identifier.
            ostream << std::setw(4) << counter;   // Residue sequence number.
            ostream << " ";
            auto r_out = p->position();
            position_t r;
            if (name == catalog->staticBP()->name() ) {
                r = r_out;
            } else {
                r = pbc.placeInside(r_out);
            }
            ostream << std::setw(8) << r[0] << std::setw(8) << r[1] << std::setw(8) << r[2];
            ostream << std::endl;
            counter += 1;
        }
        ostream << std::endl;
    });
    ostream.close();
}

int main(int argc, char *argv[]) {
    util::Logger logger{"simploce::s-pdb::main"};
    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnInputParticleSystem{"in.ps"};
    std::string fnOutputPDB{"out.pdb"};
    bool isCoarseGrained{false};

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
    if (vm.count("verbose") ) {
        util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    }

    // Read particle specifications.
    auto catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);

    // Read particle system.
    auto particleSystem = factory::particleSystem(fnInputParticleSystem, catalog, isCoarseGrained);
    logger.info("Read particle system from '" + fnInputParticleSystem + "'.");

    // Write PDB.
    toPDB(particleSystem, fnOutputPDB, catalog);
    logger.info("Wrote output PDB to: " + fnOutputPDB);

    return EXIT_SUCCESS;
}