/* 
 * File:   s-gr.cpp
 * Author: ajuffer
 *
 * Created on September 22, 2019, 3:15 PM
 */

#include "simploce/simulation/sall.hpp"
#include "simploce/simulation/sconf.hpp"
#include "simploce/analysis/gr.hpp"
#include "simploce/analysis/analysis.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/util/file.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>

namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[])
{
  using gr_t = Gr<Bead>;
  using gr_ptr_t = std::shared_ptr<gr_t>;
  using analysis_t = Analysis<Bead>;
  
  std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
  std::string fnTrajectory{"trajectory.dat"};
  std::string fnInputModel{"in.model"};
  std::string fnResults{"gr.dat"};
  std::string specName1{"Na+"};
  std::string specName2{"Cl-"};
  real_t dr{0.1};  // Bin size.

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
     "spec-name-1", po::value<std::string>(&specName1),
     "First of two particle specifications names"
    )
    (
     "spec-name-2", po::value<std::string>(&specName2),
     "Second of two particle specifications names"
    )
    (
     "bin-size", po::value<real_t>(&dr),
     "Bin size of g(r). Default is 0.1 nm."
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
    specName1 = vm["name-spec-1"].as<std::string>();
  }
  if ( vm.count("name-spec-2") ) {
    specName2 = vm["name-spec-2"].as<std::string>();
  }

  // Read particle specifications.
  spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
  std::clog << "Particle specifications: " << std::endl;
  std::clog << *catalog << std::endl;
  
  sim_model_fact_ptr_t simModelFactory = factory::simulationModelFactory(catalog);
  std::ifstream istream;
  file::open_input(istream, fnInputModel);
  cg_sim_model_ptr_t sm = simModelFactory->readCoarseGrainedFrom(istream);
  istream.close();
  std::clog << *sm << std::endl;

  box_ptr_t box = sm->box();
  bc_ptr_t bc = sm->boundaryCondition();
  gr_ptr_t gr = gr_t::create(dr, specName1, specName2, box, bc);
  
  analysis_t analysis(sm, gr);
  file::open_input(istream, fnTrajectory);
  analysis.perform(istream);
  istream.close();

  // Write results.
  auto results = gr->results();
  std::ofstream ostream;
  file::open_output(ostream, fnResults);
  ostream.setf(std::ios::scientific);
  ostream.precision(conf::PRECISION);
  for ( auto pair : results) {
    ostream << std::setw(conf::WIDTH) << pair.first;
    ostream << conf::SPACE << std::setw(conf::WIDTH) << pair.second << std::endl;
  }
  ostream.flush();
  ostream.close();
    
  std::clog << std::endl;
  std::clog << "g(r) written to '" << fnResults << "'." << std::endl;
  
  return 0;
}
