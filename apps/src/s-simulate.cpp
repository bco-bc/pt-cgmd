/*
 * Coarse-grained molecular dynamics simulation.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 10, 2019, 2:58 PM
 */

#include "simploce/conf/s-conf.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/simulation.hpp"
#include "simploce/potentials/force-field.hpp"
#include "simploce/particle/protonatable-coarse-grained.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/program-options.hpp"
#include "boost/program_options.hpp"
#include <cstdlib>
#include <string>
#include <iostream>

namespace po = boost::program_options;
using namespace simploce;
using namespace simploce::param;

int main(int argc, char *argv[]) {
    util::Logger logger("simploce::s-simulate::main");
    logger.trace("Entering.");
    util::Logger::changeLogLevel(util::Logger::LOGWARN);

    std::ofstream trajectory, data;
    try {
        std::string fnTrajectory{"trajectory.dat"};
        std::string fnSimulationData{"simulation.dat"};

        // Some default parameters
        auto param = factory::simulationParameters();

        auto nSteps = param->get<std::size_t>("simulation.nsteps");   // Number of steps.
        auto nWrite = param->get<std::size_t>("simulation.nwrite");   // Number of steps between writing
                                                                                       // to simulation data file
                                                                                       // and trajectory.
        auto nPairLists = param->get<std::size_t>("simulation.npairlists");  // Number of steps between
                                                                                              // updating pair list.
        auto timestep = param->get<real_t>("simulation.timestep");
        auto temperature = param->get<real_t>("simulation.temperature");
        auto gamma = param->get<real_t>("simulation.displacer.lvv.gamma"); // Damping rate.
        std::string displacerType = conf::VELOCITY_VERLET;                        // Default displacer.
        auto cutoffSR = param->get<real_t>("simulation.forces.cutoffSR");         // Cutoff distance for short range interactions.
        auto cutoffLR = param->get<real_t>("simulation.forces.cutoffLR");         // Cutoff distance for short range interactions.

        bool pbc_1{false};                               // Apply 1D-PBC.
        char direction{'z'};                             // Direction.
        real_t range{0.1};                               // Range for a random number sampled in MC.
        std::size_t nScaleVelocities{0};

        po::options_description usage("Usage");
        util::addStandardOptions(usage);
        usage.add_options() (
                "fn-trajectory", po::value<std::string>(&fnTrajectory),
                "Output file name trajectory. Default is 'trajectory.dat'."
        )(
                "fn-sim-data", po::value<std::string>(&fnSimulationData),
                "Output file name of simulation data. Default is 'simulation.dat'."
        )(
                "cutoff-distance-sr",
                po::value<real_t>(&cutoffSR),
                "Cutoff distance for short-ranged forces."
        )(
                "cutoff-distance-lr",
                po::value<real_t>(&cutoffLR),
                "Cutoff distance for long-ranged forces."
        )(
                "temperature,T",
                po::value<real_t>(&temperature),
                "Temperature (K). Default is 298.15 K."
        )(
                "scale-velocities",
                "Scale velocities during a warmup or initialization phase."
        )(
                "number-of-steps-between-scaling-velocities",
                po::value<std::size_t>(&nScaleVelocities),
                "Number of steps between scaling velocities"
        )(
                "damping-rate,g",
                po::value<real_t>(&gamma),
                "Damping rate (ps^-1) for Langevin/DPD heat bath."
        )(
                "number-of-steps,n",
                po::value<std::size_t>(&nSteps),
                "Number of steps. Default is 1000."
        )(
                "number-of-steps-between-save",
                po::value<std::size_t>(&nWrite),
                "Number of steps between writing to simulation data file and trajectory."
        )(
                "number-of-steps-between-pair-list-update",
                po::value<std::size_t>(&nPairLists),
                "Number of steps between updating pair list."
        )(
                "time-step,t",
                po::value<real_t>(&timestep),
                "Time step (ps). Default is 0.020 ps or 20 fs."
        )(
                "displacer,d",
                po::value<std::string>(&displacerType),
                "Displacer specification. Default is 'vv' (Velocity Verlet). "
                "Other choices are 'mc' (Monte Carlo), 'lf' (leapFrog), "
                "'lvv' (Langevin Velocity Verlet), "
                "'pt-lvv' (Langevin Velocity Verlet with Proton Transfer, NOT TESTED), "
                "'mvv-dpd' (Modified Velocity Verlet (DPD)), "
                "'s1-dpd' Splitting (DPD)"
        )(
                "range",
                po::value<real_t>(&range),
                "Length of interval from which random number is selected in MC."
        )(
                "include-external-potentials,e",
                "Include -external- potentials for force calculations."
        )(
                "pbc-1",
                po::value<char>(&direction),
                "Apply periodicity in the given direction only. Select one of 'x', y', and 'z'."
        );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << "Perform a simulation" << std::endl;
            std::cout << usage << "\n";
            return 0;
        }
        util::verbose(vm);

        // First read simulation parameters, if provided.
        util::getParameters(vm, param);
        if (vm.count("cutoff-distance-sr")) {
            cutoffSR = vm["cutoff-distance-sr"].as<real_t>();
            param->put<real_t>("simulation.forces.cutoffSR", cutoffSR);
        }
        if (vm.count("cutoff-distance-lr")) {
            cutoffSR = vm["cutoff-distance-lr"].as<real_t>();
            param->put<real_t>("simulation.forces.cutoffLR", cutoffSR);
        }
        if (vm.count("fn-trajectory")) {
            fnTrajectory = vm["fn-trajectory"].as<std::string>();
        }
        if (vm.count("fn-sim-data")) {
            fnSimulationData = vm["fn-sim-data"].as<std::string>();
        }
        if (vm.count("temperature")) {
            temperature = vm["temperature"].as<real_t>();
            param->put<real_t>("simulation.temperature", temperature);
        }
        if (vm.count("scale-velocities")) {
            param->put<bool>("simulation.scale-velocities", true);
        }
        if (vm.count("number-of-steps-between-scaling-velocities")) {
            nScaleVelocities = vm["number-of-steps-between-scaling-velocities"].as<std::size_t>();
            param->put("simulation.nscale-velocities", nScaleVelocities);
        }
        if (vm.count("damping-rate")) {
            gamma = vm["damping-rate"].as<real_t>();
            param->put<real_t>("simulation.displacer.lvv.gamma", gamma);
            param->put<real_t>("simulation.displacer.dpd.gamma", gamma);
        }
        if (vm.count("number-of-steps")) {
            nSteps = vm["number-of-steps"].as<std::size_t>();
            param->put<std::size_t>("simulation.nsteps", nSteps);
        }
        if (vm.count("number-of-steps-between-save")) {
            nWrite = vm["number-of-steps-between-save"].as<std::size_t>();
            param->put<std::size_t>("simulation.nwrite", nWrite);
        }
        if (vm.count("number-of-steps-between-pair-list-update")) {
            nPairLists = vm["number-of-steps-between-pair-list-update"].as<std::size_t>();
            param->put<std::size_t>("simulation.npairlists", nPairLists);
        }
        if (vm.count("time-step")) {
            timestep = vm["time-step"].as<real_t>();
            param->put<real_t>("simulation.timestep", timestep);
        }
        if (vm.count("displacer")) {
            displacerType = vm["displacer"].as<std::string>();
        }
        if (vm.count("range")) {
            range = vm["range"].as<real_t>();
            param->put<real_t>("simulation.displacer.mc.range", range);
        }
        if (vm.count("include-external-potentials")) {
            param->put<bool>("simulation.forces.include-external", true);
        }
        if (vm.count("pbc-1")) {
            pbc_1 = true;
            direction = vm["pbc-1"].as<char>();
        }

        logger.info("Simulation parameters:");
        std::cout << *param << std::endl;

        // Read particle system.
        auto particleSystem = util::getParticleSystem(vm, param);

       // Force field
        logger.info("Force field:");
        auto forceField = util::getForceField(vm);
        std::cout << *forceField << std::endl;

        // Interactor
        bc_ptr_t bc =  factory::pbc(particleSystem->box());
        if (pbc_1) {
            logger.info(std::to_string(direction) + ": Applying PBC in this direction only.");
            bc = factory::pbc1dBB(particleSystem->box(),
                                  Direction::valueOf(direction));
        } else {
            logger.info("Applying PBC in all directions.");
        }
        auto interactor = factory::interactor(param, forceField, bc);

        // Get the displacer.
        auto dpdUnits = factory::dpdUnits(1.0, 1.0, 1.0);
        auto displacer = factory::displacer(displacerType, param, interactor, bc, dpdUnits);

        // Simulate
        Simulation simulation(param, particleSystem, displacer, bc);
        util::open_output_file(trajectory, fnTrajectory);
        util::open_output_file(data, fnSimulationData);
        simulation.perform(trajectory, data);
        trajectory.close();
        data.close();

        // Write output particle system.
        util::writeParticleSystem(vm, particleSystem);

    } catch (std::exception &exception) {
        trajectory.flush();
        trajectory.close();
        data.flush();
        data.close();
        std::cerr << "ERROR: " << exception.what() << std::endl;
        logger.trace("Leaving.");
        return EXIT_FAILURE;
    }

    logger.trace("Leaving.");
    return EXIT_SUCCESS;
}
