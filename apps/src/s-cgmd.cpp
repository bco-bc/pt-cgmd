/*
 * Coarse-grained molecular dynamics simulation.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on September 10, 2019, 2:58 PM
 */

#include "simploce/simulation/s-types.hpp"
#include "simploce/simulation/s-conf.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/simulation/mc.hpp"
#include "simploce/simulation/simulation.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/p-factory.hpp"
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

/*
 * Coarse grained molecular dynamics.
 */
int main(int argc, char *argv[]) {
    util::Logger logger("s-cgmd::main");

    std::ofstream trajectory, data, stream;
    cg_sim_model_ptr_t model;
    std::string fnOutputModel{"out.ParticleModelSpecification"};
    try {
        std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
        std::string fnTrajectory{"trajectory.dat"};
        std::string fnSimulationData{"sim-data.dat"};
        std::string fnInputModel{};

        std::size_t n_steps = 1000;
        std::size_t n_write = 10;                         // Number of steps between writing to simulation
        // data file and trajectory.
        std::size_t npairlist = 10;                      // Number of steps between updating pair list.
        real_t timestep{0.020};                          // 0.001 ps = 1 fs. Default is 20 fs.
        real_t boxSize{6.0};                             // nm.
        real_t molarity{0.1};                            // mol/l
        real_t density{997.0479};                        // kg/m^3
        real_t temperature{298.15};                      // K.
        real_t gamma{1.0};                               // ps^-1
        std::size_t nmaxPolWaters = 1000000;             // Maximum number of polarizable waters
        // (groups).
        std::string modelType{conf::POLARIZABLE_WATER};             // Coarse grained polarizable water.
        bool doMC = false;
        // bool protonatable = true;
        std::string displacerId = conf::LANGEVIN_VELOCITY_VERLET;   // Displacer.
        real_t fc{100.0};                                // Force constant harmonic potential.
        length_t Rref{0.4};                           // Reference distance harmonic potential.
        length_t R0{0.5};                             // Initial distance between particles undergoing
        // harmonic motion.

        po::options_description usage("Usage");
        usage.add_options()
                (
                        "fn-particle-spec-catalog",
                        po::value<std::string>(&fnParticleSpecCatalog),
                        "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
                )
                (
                        "fn-input-ParticleModelSpecification", po::value<std::string>(&fnInputModel),
                        "Input file name ParticleModelSpecification. If provided, this ParticleModelSpecification will be employed for the simulation, "
                        "instead of generating a ParticleModelSpecification from scratch."
                )

                (
                        "fn-trajectory", po::value<std::string>(&fnTrajectory),
                        "Output file name trajectory. Default is 'trajectory.dat'."
                )
                (
                        "fn-output-ParticleModelSpecification", po::value<std::string>(&fnOutputModel),
                        "Output file name ParticleModelSpecification. Default is 'out.ParticleModelSpecification'."
                )
                (
                        "fn-sim-data", po::value<std::string>(&fnSimulationData),
                        "Output file name of simulation data. Default is 'sim-data.dat'."
                )

                (
                        "box-size", po::value<real_t>(&boxSize),
                        "Box (cube, unit cell) size (nm). Default is 6.0 nm."
                )
                (
                        "molarity", po::value<real_t>(&molarity),
                        "Molarity of NaCl electrolyte (mol/l). Default is 0.1 M."
                )
                (
                        "density", po::value<real_t>(&density),
                        "Water density. Default is 997.0479 kg/m^3."
                )

                (
                        "temperature", po::value<real_t>(&temperature),
                        "Temperature (K). Default is 298.15 K."
                )
                (
                        "damping-rate (Langevin heat bath)", po::value<real_t>(&gamma),
                        "Damping rate (ps^-1). Default is 1.0 ps^-1."
                )
                (
                        "number-of-steps", po::value<std::size_t>(&n_steps),
                        "Number of steps. Default is 1000."
                )
                (
                        "number-of-steps-between-save", po::value<std::size_t>(&n_write),
                        "Number of steps between writing to simulation data file and trajectory. Default is 10."
                )
                (
                        "number-of-steps-between-pairlist-update", po::value<std::size_t>(&npairlist),
                        "Number of steps between updating pair list. Default is 10."
                )
                (
                        "time-step", po::value<real_t>(&timestep),
                        "Time step (ps). Default is 0.020 ps or 20 fs."
                )

                (
                        "ParticleModelSpecification-type", po::value<std::string>(&modelType),
                        "Type of ParticleModelSpecification or system. Default is 'pol-water'. "
                        "Other choices: 'acid-base-solution', 'electrolyte', 'lj-fluid', 'hp'. "
                        "The latter creates a simulation ParticleModelSpecification that consists of a single particle group with a "
                        "single bond whose two particles undergo simple harmonic motion. This ParticleModelSpecification is "
                        "mostly used for testing purposes."
                )
                (
                        "displacer", po::value<std::string>(&displacerId),
                        "Displacer specification. Default is 'lvv' (Langevin Velocity Verlet). Other choices are"
                        "'lf' (leapFrog), 'vv' (Velocity Verlet), and "
                        "'pt-lvv' (Langevin Velocity Verlet with Proton Transfer)"
                )
                (
                        "max-number-of-water-groups", po::value<std::size_t>(&nmaxPolWaters),
                        "Maximum number of water groups. Default is 1000000."
                )
                (
                        "monte-carlo",
                        "Perform a Monte Carlo simulation."
                )
                (
                        "help", "This help message"
                );

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, usage), vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << "Coarse-grained molecular dynamics simulation" << std::endl;
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
        if (vm.count("fn-sim-data")) {
            fnSimulationData = vm["fn-sim-data"].as<std::string>();
        }
        if (vm.count("fn-output-ParticleModelSpecification")) {
            fnOutputModel = vm["fn-output-ParticleModelSpecification"].as<std::string>();
        }
        if (vm.count("box-size")) {
            boxSize = vm["box-size"].as<real_t>();
        }
        if (vm.count("molarity")) {
            molarity = vm["molarity"].as<real_t>();
        }
        if (vm.count("density")) {
            density = vm["density"].as<real_t>();
        }
        if (vm.count("temperature")) {
            temperature = vm["temperature"].as<real_t>();
        }
        if (vm.count("damping-rate")) {
            gamma = vm["damping-rate"].as<real_t>();
        }
        if (vm.count("number-of-steps")) {
            n_steps = vm["number-of-steps"].as<std::size_t>();
        }
        if (vm.count("number-of-steps-between-save")) {
            n_write = vm["number-of-steps-between-save"].as<std::size_t>();
        }
        if (vm.count("number-of-steps-between-pairlist-update")) {
            npairlist = vm["number-of-steps-between-pairlist-update"].as<std::size_t>();
        }
        if (vm.count("time-step")) {
            timestep = vm["time-step"].as<real_t>();
        }
        if (vm.count("ParticleModelSpecification-type")) {
            modelType = vm["ParticleModelSpecification-type"].as<std::string>();
        }
        if (vm.count("displacer")) {
            displacerId = vm["displacer"].as<std::string>();
        }
        if (vm.count("max-number-of-water-groups")) {
            nmaxPolWaters = vm["max-number-of-water-groups"].as<std::size_t>();
        }
        if (vm.count("monte-carlo")) {
            doMC = true;
        }

        // Simulation parameters
        sim_param_t param;
        param.add<std::size_t>("n_steps", n_steps);
        param.add<std::size_t>("n_write", n_write);
        param.add<real_t>("temperature", temperature);
        param.add<real_t>("timestep", timestep);
        param.add<real_t>("gamma", gamma);
        param.add<std::size_t>("npairlists", npairlist);
        logger.info("Simulation parameters:");
        std::cout << param << std::endl;

        // Read particle specifications.
        spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
        std::clog << "Particle specifications: " << std::endl;
        std::clog << *catalog << std::endl;

        // Read or set up new simulation ParticleModelSpecification.
        sim_model_fact_ptr_t simModelFactory = factory::simulationModelFactory(catalog);
        if (fnInputModel.empty()) {
            box_ptr_t box = std::make_shared<box_t>(boxSize);
            if (modelType == conf::POLARIZABLE_WATER) {
                model = simModelFactory->polarizableWater(box, density, temperature, nmaxPolWaters);
            } else if (modelType == conf::ACID_BASE_SOLUTION) {
             /* ParticleModelSpecification = simModelFactory->formicAcidSolution(box,
                                                            density,
                                                            molarity,
                                                            temperature,
                                                            protonatable);*/
            } else if (modelType == conf::ELECTROLYTE) {
                model = simModelFactory->electrolyte(box, molarity, temperature);
            } else if (modelType == conf::LJ_FLUID) {
                model = simModelFactory->ljFluid(box, density, temperature);
            } else if (modelType == conf::HP) {
                model = simModelFactory->harmonic(R0, Rref, fc);
            } else {
                throw std::domain_error(modelType + ": No such simulation ParticleModelSpecification type available (yet).");
            }
            factory::changeDisplacer(displacerId, model);
            std::clog << std::endl;
            std::clog << "Created molecular ParticleModelSpecification:" << std::endl;

            std::ofstream str;
            util::open_output_file(str, "created.ParticleModelSpecification");
            str << *model << std::endl;
            str.close();

            std::clog << "Created ParticleModelSpecification saved in 'created.ParticleModelSpecification'." << std::endl;
            std::clog << std::endl;
        } else {
            std::ifstream str;
            util::open_input_file(str, fnInputModel);
            model = simModelFactory->readCoarseGrainedFrom(str);
            stream.close();
            std::clog << "Read ParticleModelSpecification from input file '" << fnInputModel << "'." << std::endl;
        }

        // Simulate.
        util::open_output_file(trajectory, fnTrajectory);
        util::open_output_file(data, fnSimulationData);
        if (doMC) {
            MC<Bead> mc(model);
            mc.perform(param, trajectory, data);
        } else {
            Simulation<Bead> simulation(model);
            simulation.perform(param, trajectory, data);
        }
        trajectory.close();
        data.close();

        // Write output ParticleModelSpecification.
        util::open_output_file(stream, fnOutputModel);
        stream << *model << std::endl;
        stream.close();
        std::clog << "Wrote ParticleModelSpecification to output file '" << fnOutputModel << "'." << std::endl;

        return 0;
    } catch (std::exception &exception) {
        trajectory.flush();
        trajectory.close();
        data.flush();
        data.close();
        if (model) {
            util::open_output_file(stream, fnOutputModel);
            stream << *model << std::endl;
            stream.close();
        }
        std::cerr << "ERROR: " << exception.what() << std::endl;
    }
}
