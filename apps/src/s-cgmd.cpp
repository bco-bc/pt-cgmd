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
    std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
    std::string fnTrajectory{"trajectory.dat"};
    std::string fnSimulationData{"sim-data.dat"};
    std::string fnInputModel{};
    std::string fnOutputModel{"out.model"};

    std::size_t nsteps = 1000;
    std::size_t nwrite = 10;
    real_t timestep{0.020};                          // 0.001 ps = 1 fs. Default is 20 fs.
    real_t boxSize{6.0};                             // nm.
    real_t molarity{0.1};                            // mol/l
    real_t density{997.0479};                        // kg/m^3
    real_t temperature{298.15};                      // K.
    real_t gamma{0.50};                              // ps^-1
    std::string modelType{conf::POLARIZABLE_WATER};  // Coarse grained polarizable water.

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
       "Damping rate (ps^-1). Default is 0.5 ps^-1."
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
       "Other choices: acid-base-solution'.")
      
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
    cg_sim_model_ptr_t model;
    if ( fnInputModel.empty() ) {
      box_ptr_t box = std::make_shared<box_t>(boxSize);      
      if ( modelType == conf::POLARIZABLE_WATER ) {
	model = simModelFactory->polarizableWater(box, density, temperature);
      } else if ( modelType == conf::ACID_BASE_SOLUTION ) {
	model = simModelFactory->formicAcidSolution(box, density, molarity, temperature, true);
      } else {
	throw std::domain_error(modelType + ": No such simulation model type available (yet).");
      }
      std::clog << std::endl;
      std::clog << "Created molecular model:" << std::endl;
      
      std::cout << *model << std::endl;
    } else {
      std::ifstream stream;
      file::open_input(stream, fnInputModel);
      model = simModelFactory->readCoarseGrainedFrom(stream);
      std::clog << "Read model from input file '" << fnInputModel << "'." << std::endl;
      stream.close();
    }

    // Simulate.
    Simulation<Bead> simulation(model);
    std::ofstream traj, data;
    file::open_output(traj, fnTrajectory);
    file::open_output(data, fnSimulationData);
    simulation.perform(param, traj, data);
    traj.close();
    data.close();
    
    // Write output model.
    std::ofstream stream;
    file::open_output(stream, fnOutputModel);
    stream << *model << std::endl;
    stream.close();
    
    return 0;    
}
