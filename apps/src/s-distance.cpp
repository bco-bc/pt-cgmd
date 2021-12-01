/*
 * Author: André H. Juffer.
 * Created on 28/11/2021, 14:56.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/analysis/distance.hpp"
#include "simploce/analysis/analysis.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/logger.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <cstdlib>

namespace po = boost::program_options;
using namespace simploce;
using namespace simploce::param;

int main(int argc, char *argv[]) {
    util::Logger logger("simploce::s-distances::main");

    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnTrajectory{"trajectory.dat"};
    std::string fnInputParticleSystem{"in.ps"};
    std::string fnResults{"distances.dat"};
    std::string idi{"1"};
    std::string idj{"2"};
    real_t dr{0.01};              // Bin size, nm.
    std::size_t nSkip = 0;        // Number of states in trajectory to skip
                                  // before executing analysis.
    bool isCoarseGrained{false};

    po::options_description usage("Usage");
    usage.add_options() (
        "fn-particle-spec-catalog,s",
        po::value<std::string>(&fnParticleSpecCatalog),
        "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
    )(
        "fn-input-particle-system,i", po::value<std::string>(&fnInputParticleSystem),
        "Input file name particle system. Default is 'in.ps'."
    )(
        "fn-trajectory",
        po::value<std::string>(&fnTrajectory),
        "Input file name trajectory. Default is 'trajectory.dat'."
    )(
        "id-1,1",
        po::value<std::string>(&idi),
        "First of two particle identifiers"
    )(
        "id-2,2", po::value<std::string>(&idj),
        "Second of two particle identifiers"
    )(
        "coarse-grained,c",
        "Input is a coarse-grained description."
    )(
        "skip-number-of-states,k",
        po::value<std::size_t>(&nSkip),
        "Number of states in trajectory to skip before executing analysis"
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
        std::cout << "Distance between two particles" << std::endl;
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

    if (vm.count("fn-trajectory")) {
        fnTrajectory = vm["fn-trajectory"].as<std::string>();
    }

    if (vm.count("id-1")) {
        idi = vm["id-1"].as<std::string>();
    }
    if (vm.count("id-2")) {
        idj = vm["id-2"].as<std::string>();
    }
    if (vm.count("skip-number-of-states")) {
        nSkip = vm["skip-number-of-states"].as<std::size_t>();
    }
    if (vm.count("verbose") ) {
        util::Logger::changeLogLevel(util::Logger::LOGDEBUG);
    }

    // Analysis parameters
    a_param_ptr_t analysisParameters = factory::simulationParameters();
    analysisParameters->add<std::size_t>("analysis.trajectory.nskip", nSkip);
    analysisParameters->add<std::string>("analysis.distance.id_i", idi);
    analysisParameters->add<std::string>("analysis.distance.id_j", idj);
    logger.info("Analysis parameters:");
    std::cout << *analysisParameters << std::endl;

    // Read particle specifications.
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);

    // Read particle system.
    auto particleSystem = factory::particleSystem(fnInputParticleSystem, catalog, isCoarseGrained);
    logger.info("Read particle system from '" + fnInputParticleSystem + "'.");

    // Analyzer.
    bc_ptr_t bc = factory::boundaryCondition(particleSystem->box());
    d_ptr_t distances = Distance::create(analysisParameters, bc);

    // Perform analysis.
    Analysis analysis(particleSystem, analysisParameters, distances);
    std::ifstream trajectory;
    util::open_input_file(trajectory, fnTrajectory);
    analysis.perform(trajectory);
    trajectory.close();

    // Write results.
    auto results = distances->results();
    std::ofstream ostream;
    util::open_output_file(ostream, fnResults);
    ostream.setf(std::ios::scientific);
    ostream.precision(conf::PRECISION);
    for (auto pair: results) {
        ostream << std::setw(conf::REAL_WIDTH) << pair.first;
        ostream << conf::SPACE << std::setw(conf::REAL_WIDTH) << pair.second << std::endl;
    }
    ostream.flush();
    ostream.close();
    logger.info("Distances written to '" + fnResults + "'.");

    // Done
    return EXIT_SUCCESS;
}

