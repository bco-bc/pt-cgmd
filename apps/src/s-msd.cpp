/*
 * Diffusion constant.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on November 15, 2019, 1:57 PM
 */

#include "simploce/analysis/diffusion.hpp"
#include "simploce/analysis/analysis.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/units/units-mu.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>

namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[]) {
    using analyzer_t = Diffusion<Bead>;
    using analyzer_ptr_t = std::shared_ptr<analyzer_t>;
    using analysis_t = Analysis<Bead>;

    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnTrajectory{"trajectory.dat"};
    std::string fnInputModel{"in.ParticleModelSpecification"};
    std::string fnResults{"msd.dat"};
    real_t dt{10 * 0.02};         // Time interval between successive states in trajectory.
    real_t t0{0.0};               // Start time.
    real_t tau(10.0);             // Window time interval for mean square deviation.
    std::size_t nSkip = 0;        // Number of states in trajectory to skip
                                  // before executing analysis.

    po::options_description usage("Usage");
    usage.add_options()
            (
                    "fn-particle-spec-catalog",
                    po::value<std::string>(&fnParticleSpecCatalog),
                    "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
            )
            (
                    "fn-ParticleModelSpecification", po::value<std::string>(&fnInputModel),
                    "Input file name ParticleModelSpecification. Default is 'in.ParticleModelSpecification'.")
            (
                    "fn-trajectory", po::value<std::string>(&fnTrajectory),
                    "Input file name trajectory. Default is 'trajectory.dat'."
            )
            (
                    "time-interval", po::value<real_t>(&dt),
                    "Time interval between successive states in trajectory."
                    "Default is 0.2 ps."
            )
            (
                    "start-time", po::value<real_t>(&t0),
                    "Start time. Default is 0. If a number of states in the trajectory must be skipped, "
                    "then the start time would change accordingly."
            )
            (
                    "skip-number-of-states", po::value<std::size_t>(&nSkip),
                    "Number of states in trajectory to skip before executing analysis"
            )
            (
                    "help", "Help message"
            );

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, usage), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << std::endl;
        std::cout << usage << "\n";
        return 0;
    }
    if (vm.count("fn-particle-spec-catalog")) {
        fnParticleSpecCatalog = vm["fn-particle-spec-catalog"].as<std::string>();
    }
    if (vm.count("fn-input-ParticleModelSpecification")) {
        fnInputModel = vm["fn-input-ParticleModelSpecification"].as<std::string>();
    }

    if (vm.count("fn-trajectory")) {
        fnTrajectory = vm["fn-trajectory"].as<std::string>();
    }

    if (vm.count("time-interval")) {
        dt = vm["time-interval"].as<real_t>();
    }
    if (vm.count("skip-number-of-states")) {
        nSkip = vm["skip-number-of-states"].as<std::size_t>();
    }

    // Simulation parameters
    sim_param_t param;
    param.add<std::size_t>("nSkip", nSkip);

    // Read particle specifications.
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
    std::clog << "Particle specifications: " << std::endl;
    std::clog << *catalog << std::endl;

    sim_model_fact_ptr_t simModelFactory = factory::simulationModelFactory(catalog);
    std::ifstream istream;
    util::open_input_file(istream, fnInputModel);
    cg_sim_model_ptr_t sm = simModelFactory->readCoarseGrainedFrom(istream);
    istream.close();
    std::clog << "Number of particles: " << sm->size() << std::endl;

    // Perform analysis.
    analyzer_ptr_t analyzer = analyzer_t::create(dt, tau);
    analysis_t analysis(sm, analyzer);
    util::open_input_file(istream, fnTrajectory);
    analysis.perform(param, istream);
    istream.close();

    // Write results.
    auto results = analyzer->results();
    std::ofstream ostream;
    util::open_output_file(ostream, fnResults);
    ostream.setf(std::ios::scientific);
    ostream.precision(conf::PRECISION);
    for (auto &result: results) {
        ostream << std::setw(conf::REAL_WIDTH) << result.first
                << std::setw(conf::REAL_WIDTH) << result.second
                << std::endl;
    }
    ostream.flush();
    ostream.close();

    std::clog << "MSD was written to output file " << fnResults << "'." << std::endl;

    return 0;
}
