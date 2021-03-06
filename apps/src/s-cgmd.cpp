/* 
 * File:   s-cgmd.cpp
 * Author: ajuffer
 *
 * Created on September 10, 2019, 2:58 PM
 */

#include "simploce/simulation/sall.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
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
int main(int argc, char* argv[]) 
{
  std::ofstream traj, data, stream;
  cg_sim_model_ptr_t model;
  std::string fnOutputModel{"out.model"};
  try {
    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnTrajectory{"trajectory.dat"};
    std::string fnSimulationData{"sim-data.dat"};
    std::string fnInputModel{};

    std::size_t nsteps = 1000;
    std::size_t nwrite = 10;
    real_t timestep{0.020};                          // 0.001 ps = 1 fs. Default is 20 fs.
    real_t boxSize{6.0};                             // nm.
    real_t molarity{0.1};                            // mol/l
    real_t density{997.0479};                        // kg/m^3
    real_t temperature{298.15};                      // K.
    real_t gamma{1.0};                               // ps^-1
    std::size_t nmaxPolWaters = 1000000;             // Maximum number of polarizable waters
                                                     // (groups).
    std::string modelType{conf::POLARIZABLE_WATER};  // Coarse grained polarizable water.
    bool mc = false;
    bool protonatable = true;
    std::string displacerId =
      conf::LANGEVIN_VELOCITY_VERLET;                // Displacer.
    real_t fc{100.0};                                // Force constant harmonic potential.
    length_t Rref{0.4};                              // Reference distance harmonic potential.
    length_t R0{0.5};                                // Initial distance between particles undergoing
                                                     // harmonic motion.

    po::options_description usage("Usage");
    usage.add_options()
      (
       "fn-particle-spec-catalog",
       po::value<std::string>(&fnParticleSpecCatalog),
       "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
       )
      (
       "fn-input-model",  po::value<std::string>(&fnInputModel),
       "Input file name model. If provided, this model will be employed for the simulation, "
       "instead of generating a model from scratch."
       )
      
      (
       "fn-trajectory",  po::value<std::string>(&fnTrajectory),
       "Output file name trajectory. Default is 'trajectory.dat'."
       )
      (
       "fn-output-model",  po::value<std::string>(&fnOutputModel),
       "Output file name model. Default is 'out.model'."
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
       "number-of-steps", po::value<std::size_t>(&nsteps),
       "Number of steps. Default is 1000."
      )
      (
       "number-steps-between-save", po::value<std::size_t>(&nwrite),
       "Number of steps between saving state. Default is 10."
      )
      (
       "time-step-size", po::value<real_t>(&timestep),
       "Time step size (ps). Default is 0.020 ps or 20 fs."
       )
      
      (
       "model-type", po::value<std::string>(&modelType),
       "Type of model or system. Default is 'pol-water'. "
       "Other choices: 'acid-base-solution', 'electrolyte', 'lj-fluid', 'hp'. "
       "The latter creates a simulation model that consists of a single particle group with a "
       "single bond whose two particles undergo simple harmonic motion. This model is "
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
       "help", "Help message"
      )
      ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, usage), vm);
    po::notify(vm);
    
    if ( vm.count("help") || argc == 1) {
        std::cout << std::endl;
        std::cout << usage << "\n";
        return 0;
    }

    if ( vm.count("fn-particle-spec-catalog") ) {
      fnParticleSpecCatalog = vm["fn-particle-spec-catalog"].as<std::string>();
    }
    if ( vm.count("fn-input-model") ) {
      fnInputModel = vm["fn-input-model"].as<std::string>();
    }
    if ( vm.count("fn-trajectory") ) {
      fnTrajectory = vm["fn-trajectory"].as<std::string>();
    }
    if ( vm.count("fn-sim-data") ) {
      fnSimulationData = vm["fn-sim-data"].as<std::string>();
    }
    if ( vm.count("fn-output-model") ) {
      fnOutputModel = vm["fn-output-model"].as<std::string>();
    }
    if ( vm.count("box-size") ) {
      boxSize =  vm["box-size"].as<real_t>();
    }
    if ( vm.count("molarity") ) {
      molarity = vm["molarity"].as<real_t>();
    }
    if ( vm.count("density") ) {
      density = vm["density"].as<real_t>();
    }
    if ( vm.count("temperature") ) {
      temperature = vm["temperature"].as<real_t>();
    }
    if ( vm.count("damping-rate") ) {
      gamma = vm["damping-rate"].as<real_t>();
    }
    if ( vm.count("number-of-steps") ) {
      nsteps = vm["number-of-steps"].as<std::size_t>();
    }
    if ( vm.count("number-of-steps-between-save") ) {
      nwrite = vm["number-of-steps"].as<std::size_t>();
    }
    if (vm.count("time-step-size") ) {
      timestep = vm["time-step-size"].as<real_t>();
    }
    if ( vm.count("model-type") ) {
      modelType = vm["model-type"].as<std::string>();
    }
    if ( vm.count("displacer") ) {
      displacerId = vm["displacer"].as<std::string>();
    }
    if ( vm.count("max-number-of-water-groups") ) {
      nmaxPolWaters = vm["max-number-of-water-groups"].as<std::size_t>();
    }
    if ( vm.count("monte-carlo") ) {
      mc = true;
    }
    
    // Simulation parameters
    sim_param_t param;    
    param.add<std::size_t>("nsteps", nsteps);
    param.add<std::size_t>("nwrite", nwrite);
    param.add<real_t>("temperature", temperature);
    param.add<real_t>("timestep", timestep);
    param.add<real_t>("gamma", gamma);
    param.add<std::size_t>("npairlists", 10);
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << param << std::endl;
    
    // Read particle specifications.
    spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
    std::clog << "Particle specifications: " << std::endl;
    std::clog << *catalog << std::endl;

    // Read or set up new simulation model.
    sim_model_fact_ptr_t simModelFactory = factory::simulationModelFactory(catalog);
    if ( fnInputModel.empty() ) {
      box_ptr_t box = std::make_shared<box_t>(boxSize);      
      if ( modelType == conf::POLARIZABLE_WATER ) {
	model = simModelFactory->polarizableWater(box, density, temperature, nmaxPolWaters);
      } else if ( modelType == conf::ACID_BASE_SOLUTION ) {
	model = simModelFactory->formicAcidSolution(box,
						    density,
						    molarity,
						    temperature,
						    protonatable);
      } else if ( modelType == conf::ELECTROLYTE ) {
	model = simModelFactory->electrolyte(box, molarity, temperature);	
      } else if ( modelType == conf::LJ_FLUID) {
	model = simModelFactory->ljFluid(box, density, temperature);
      } else if ( modelType == conf::HP) {
	model = simModelFactory->harmonic(R0, Rref, fc);
      } else {
	throw std::domain_error(modelType + ": No such simulation model type available (yet).");
      }
      factory::changeDisplacer(displacerId, model);
      std::clog << std::endl;
      std::clog << "Created molecular model:" << std::endl;
      
      std::ofstream stream;
      file::open_output(stream, "created.model");
      stream << *model << std::endl;
      stream.close();
      
      std::clog << "Created model saved in 'created.model'." << std::endl;
      std::clog << std::endl;
    } else {
      std::ifstream stream;
      file::open_input(stream, fnInputModel);
      model = simModelFactory->readCoarseGrainedFrom(stream);
      stream.close();
      std::clog << "Read model from input file '" << fnInputModel << "'." << std::endl;
    }

    // Simulate.
    file::open_output(traj, fnTrajectory);
    file::open_output(data, fnSimulationData);
    if ( mc ) {
      MC<Bead> mc(model);
      mc.perform(param, traj, data);
    } else {
      Simulation<Bead> simulation(model);
      simulation.perform(param, traj, data);
    }
    traj.close();
    data.close();
    
    // Write output model.
    file::open_output(stream, fnOutputModel);
    stream << *model << std::endl;
    stream.close();
    std::clog << "Wrote model to output file '" << fnOutputModel << "'." << std::endl;
    
    return 0;
  } catch (std::exception& exception) {
    traj.flush();
    traj.close();
    data.flush();
    data.close();
    if ( model ) {
      file::open_output(stream, fnOutputModel);
      stream << *model << std::endl;
      stream.close();
    }
    std::cerr << exception.what() << std::endl;
  }
}
