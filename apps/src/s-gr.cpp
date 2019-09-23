/* 
 * File:   s-gr.cpp
 * Author: ajuffer
 *
 * Created on September 22, 2019, 3:15 PM
 */

#include "simploce/simulation/stypes.hpp"
#include "boost/program_options.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[])
{
  std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
  std::string fnTrajectory{"trajectory.dat"};
  std::string fnInputModel{"in.model"};
  std::string fnResults{"gr.dat"};
  std::string nameSpec1{"Na+"};
  std::string nameSpec2{"Cl-"};

  po::options_description usage("Usage");
  usage.add_options()
    (
       "fn-particle-spec-catalog",
       po::value<std::string>(&fnParticleSpecCatalog),
       "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
    )
    (
     "fn-model",  po::value<std::string>(&fnInputModel),
     "Input file name model. Default is 'in.model'.")
    (
     "fn-trajectory",  po::value<std::string>(&fnTrajectory),
     "Input file name trajectory. Default is 'trajectory.dat'."
    )
    (
     "name-spec-1", po::value<std::string>(&nameSpec1), "First of two particle specifications names"
    )
    (
     "name-spec-2", po::value<std::string>(&nameSpec2), "Second of two particle specifications names"
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

  if ( vm.count("name-spec-1") ) {
    nameSpec1 = vm["name-spec-1"].as<std::string>();
  }
  if ( vm.count("name-spec-2") ) {
    nameSpec2 = vm["name-spec-2"].as<std::string>();
  }

  /*
  Factory& factory = Factory::create(fnParticleCatalog);
  particle_catalog_ptr_t particleCatalog = factory.particleCatalog();
  force_field_fact_ptr_t forceFieldFactory = factory.forceFieldFactory();
  std::ifstream istream;
  openInputFile(istream, fnModel);
  mol_model_ptr_t model = std::make_shared<MolecularModel>();
  model->readFromStream(istream, particleCatalog, forceFieldFactory);
  std::clog << "Read molecular system from '" << fnModel << "'." << std::endl;
  std::clog << *model << std::endl;

  // Set up analysis.
  particle_spec_ptr_t spec1 = particleCatalog->lookup(name1);
  particle_spec_ptr_t spec2 = particleCatalog->lookup(name2);
  store_ptr_t traj = FileStore::createForParsing(fnTrajectory);
  //analyzer_ptr_t gr = GR::make(spec1, spec2, traj);
  GR gr{spec1, spec2, traj};

  // Do analysis.
  std::clog << "Calculating g(r)..." << std::endl;
  gr.execute(model);
  std::clog << "Done." << std::endl;

  // Write results.
  auto result = gr.getResult();
  std::ofstream ostream;
  openOutputFile(ostream, fnResults);
  ostream.setf(std::ios::scientific);
  ostream.precision(PRECISION);
  for ( auto pair : result) {
    ostream << std::setw(WIDTH) << pair.first;
    ostream << SPACE << std::setw(WIDTH) << pair.second << std::endl;
  }
  ostream.flush();
  ostream.close();
    
  std::clog << std::endl;
  std::clog << "g(r) written to '" << fnResults << "'." << std::endl;
  */
  
  return 0;
}
