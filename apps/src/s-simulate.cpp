/*
 * Coarse-grained molecular dynamics simulation.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 10, 2019, 2:58 PM
 */

#include "simploce/simulation/s-types.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/simulation.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/logger.hpp"
#include <boost/program_options.hpp>
#include <cstdlib>
#include <string>
#include <iostream>

namespace po = boost::program_options;
using namespace simploce;
using namespace simploce::param;

int main(int argc, char *argv[]) {
    util::Logger logger("simploce::s-simulate::main");

    std::ofstream trajectory, data;
    try {
        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnTrajectory{"trajectory.dat"};
        std::string fnSimulationData{"simulation.dat"};
        std::string fnInputParticleSystem{"in.ps"};
        std::string fnOutputParticleSystem{"out.ps"};
        std::string fnForceField{"interaction-parameters.dat"};
        bool isCoarseGrained{false};
        bool includeExternalPotentials{false};

        std::size_t nSteps = 1000;                       // Number of steps.
        std::size_t nWrite = 10;                         // Number of steps between writing to simulation
                                                         // data file and trajectory.
        std::size_t nPairLists = 10;                     // Number of steps between updating pair list.
        real_t timestep{0.020};                          // 0.001 ps = 1 fs. Default is 20 fs.
        real_t temperature{298.15};                      // in K.
        real_t gamma{1.0};                               // in ps^-1
        std::string displacerType =
                conf::LANGEVIN_VELOCITY_VERLET;          // Displacer.

        po::options_description usage("Usage");
        usage.add_options() (
            "fn-particle-spec-catalog,s",
            po::value<std::string>(&fnParticleSpecCatalog),
            "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
        )(
            "fn-input-particle-system,i",
            po::value<std::string>(&fnInputParticleSystem),
            "Input file name particle system. Default is 'in.ps'"
        )(
            "fn-force-field,f",
            po::value<std::string>(&fnForceField),
            "Input file name interaction parameters. Default is 'interaction-parameters.dat'."
        )(
            "fn-trajectory", po::value<std::string>(&fnTrajectory),
            "Output file name trajectory. Default is 'trajectory.dat'."
        )(
            "fn-output-particle-system,o", po::value<std::string>(&fnOutputParticleSystem),
            "Output file name particle system. Default is 'out.ps'."
        )(
            "fn-sim-data", po::value<std::string>(&fnSimulationData),
            "Output file name of simulation data. Default is 'simulation.dat'."
        )(
                "coarse-grained,c",
                "Input is a coarse-grained description."
        )(
            "temperature,T", po::value<real_t>(&temperature),
            "Temperature (K). Default is 298.15 K."
        )(
            "damping-rate,g", po::value<real_t>(&gamma),
            "Damping rate (ps^-1) for Langevin heat bath. Default is 1.0 ps^-1."
        )(
            "number-of-steps,n", po::value<std::size_t>(&nSteps),
            "Number of steps. Default is 1000."
        )(
            "number-of-steps-between-save", po::value<std::size_t>(&nWrite),
            "Number of steps between writing to simulation data file and trajectory. Default is 10."
        )(
            "number-of-steps-between-pairlist-update", po::value<std::size_t>(&nPairLists),
            "Number of steps between updating pair list. Default is 10."
        )(
            "time-step,t", po::value<real_t>(&timestep),
            "Time step (ps). Default is 0.020 ps or 20 fs."
        )(
            "displacer,d", po::value<std::string>(&displacerType),
            "Displacer specification. Default is 'lvv' (Langevin Velocity Verlet). "
            "Other choices are 'mc' (Monte Carlo), 'lf' (leapFrog), 'vv' (Velocity Verlet), and "
            "'pt-lvv' (Langevin Velocity Verlet with Proton Transfer)"
        )(
            "include-external-potentials,e",
            "Include external potentials for force calculations."
        )(
            "verbose,v",
            "Verbose"
        )(
            "help,h",
            "This help message"
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << "Perform a simulation" << std::endl;
            std::cout << usage << "\n";
            return 0;
        }

        if (vm.count("fn-particle-spec-catalog")) {
            fnParticleSpecCatalog = vm["fn-particle-spec-catalog"].as<std::string>();
        }
        if (vm.count("fn-input-particle-system")) {
            fnInputParticleSystem = vm["fn-input-particle-system"].as<std::string>();
        }
        if (vm.count("fn-force-field")) {
            fnForceField = vm["fn-force-field"].as<std::string>();
        }
        if (vm.count("fn-trajectory")) {
            fnTrajectory = vm["fn-trajectory"].as<std::string>();
        }
        if (vm.count("fn-sim-data")) {
            fnSimulationData = vm["fn-sim-data"].as<std::string>();
        }
        if (vm.count("fn-output-particle-system")) {
            fnOutputParticleSystem = vm["fn-output-particle-system"].as<std::string>();
        }
        if (vm.count("coarse-grained") ) {
            isCoarseGrained = true;
        }
        if (vm.count("temperature")) {
            temperature = vm["temperature"].as<real_t>();
        }
        if (vm.count("damping-rate")) {
            gamma = vm["damping-rate"].as<real_t>();
        }
        if (vm.count("number-of-steps")) {
            nSteps = vm["number-of-steps"].as<std::size_t>();
        }
        if (vm.count("number-of-steps-between-save")) {
            nWrite = vm["number-of-steps-between-save"].as<std::size_t>();
        }
        if (vm.count("number-of-steps-between-pairlist-update")) {
            nPairLists = vm["number-of-steps-between-pairlist-update"].as<std::size_t>();
        }
        if (vm.count("time-step")) {
            timestep = vm["time-step"].as<real_t>();
        }
        if (vm.count("displacer")) {
            displacerType = vm["displacer"].as<std::string>();
        }
        if (vm.count("include-external-potentials")) {
            includeExternalPotentials = true;
        }
        if (vm.count("verbose") ) {
            logger.changeLogLevel(util::Logger::LOGTRACE);
        }

        // Simulation parameters
        sim_param_ptr_t simulationParameters = factory::simulationParameters();
        simulationParameters->put<std::size_t>("simulation.nsteps", nSteps);
        simulationParameters->put<std::size_t>("simulation.nwrite", nWrite);
        simulationParameters->put<real_t>("simulation.temperature", temperature);
        simulationParameters->put<real_t>("simulation.timestep", timestep);
        simulationParameters->put<real_t>("simulation.gamma", gamma);
        simulationParameters->put<std::size_t>("simulation.npairlists", nPairLists);
        simulationParameters->put<bool>("simulation.include-external", includeExternalPotentials);
        logger.info("Simulation parameters:");
        std::cout << *simulationParameters << std::endl;

        // Read particle specifications.
        spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);

        // Read particle system.
        auto particleSystem = factory::particleSystem(fnInputParticleSystem, catalog, isCoarseGrained);
        logger.info("Read particle system from '" + fnInputParticleSystem + "'.");

        // Force field
        auto forceField = factory::forceField(fnForceField, catalog);

        // Interactor
        auto bc = factory::boundaryCondition(particleSystem->box());
        auto interactor = factory::interactor(simulationParameters, forceField, bc);

        // Get the displacer.
        auto displacer = factory::displacer(displacerType, simulationParameters, interactor);

        // Simulate
        Simulation simulation(simulationParameters, particleSystem, displacer);
        logger.info("Simulation ongoing...");
        util::open_output_file(trajectory, fnTrajectory);
        util::open_output_file(data, fnSimulationData);
        simulation.perform(trajectory, data);
        trajectory.close();
        data.close();
        logger.info("Done.");

        // Write output particle system.
        std::ofstream stream;
        util::open_output_file(stream, fnOutputParticleSystem);
        stream << *particleSystem << std::endl;
        stream.close();
        logger.info("Output particle system written to '" + fnOutputParticleSystem + "'.");

        // Done.
        return 0;
    } catch (std::exception &exception) {
        trajectory.flush();
        trajectory.close();
        data.flush();
        data.close();
        std::cerr << "ERROR: " << exception.what() << std::endl;
    }
}
