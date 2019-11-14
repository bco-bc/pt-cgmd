/* 
 * File:   s-dipole-moment.cpp
 * Author: ajuffer
 *
 * Created on October 24, 2019, 3:15 PM
 */

#include "simploce/analysis/dipole-moment.hpp"
#include "simploce/analysis/analysis.hpp"
#include "simploce/simulation/sfactory.hpp"
#include "simploce/simulation/sim-model.hpp"
#include "simploce/simulation/sim-model-factory.hpp"
#include "simploce/simulation/sim-util.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/mu-units.hpp"
#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <utility>
#include <tuple>

namespace po = boost::program_options;
using namespace simploce;

int main(int argc, char *argv[])
{
  using analyzer_t = DipoleMoment<Bead>;
  using analyzer_ptr_t = std::shared_ptr<analyzer_t>;
  using analysis_t = Analysis<Bead>;
  
  std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
  std::string fnTrajectory{"trajectory.dat"};
  std::string fnInputModel{"in.model"};
  std::string fnResultsM{"dipole-moment.dat"};
  std::string fnResultsFM{"f-dipole-moment.dat"};
  real_t dt{10 * 0.02};         // Time interval between successive states in trajectory.
  real_t t0{0.0};               // Start time.
  std::size_t nskip = 0;        // Number of states in trajectory to skip
                                // before executing analysis.
  real_t temperature(298.15);
  real_t dm = 0.001;
  real_t mmax = 1.0;

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
     "skip-number-of-states", po::value<std::size_t>(&nskip),
     "Number of states in trajectory to skip before executing analysis"
    )
    (
     "temperature", po::value<real_t>(&temperature),
     "Temperature. Default is 298.15 K."
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

  if ( vm.count("time-interval") ) {
      dt = vm["time-interval"].as<real_t>();
  }
  if ( vm.count("skip-number-of-states") ) {
    nskip = vm["skip-number-of-states"].as<std::size_t>();
  }
  if ( vm.count("temperature") ) {
    temperature = vm["temperature"].as<real_t>();
  }

  // Simulation parameters
  sim_param_t param;
  param.add<std::size_t>("nskip", nskip);

  // Read particle specifications.
  spec_catalog_ptr_t catalog = factory::particleSpecCatalog(fnParticleSpecCatalog);
  std::clog << "Particle specifications: " << std::endl;
  std::clog << *catalog << std::endl;
  
  sim_model_fact_ptr_t simModelFactory = factory::simulationModelFactory(catalog);
  std::ifstream istream;
  file::open_input(istream, fnInputModel);
  cg_sim_model_ptr_t sm = simModelFactory->readCoarseGrainedFrom(istream);
  istream.close();
  std::clog << "Number of particles: " << sm->size() << std::endl;

  // Perform analysis.
  analyzer_ptr_t analyzer = analyzer_t::create(dt, dm, mmax, t0);
  analysis_t analysis(sm, analyzer);
  file::open_input(istream, fnTrajectory);
  analysis.perform(param, istream);
  istream.close();

  // Write results.
  auto results = analyzer->results();
  std::ofstream ostream;
  file::open_output(ostream, fnResultsM);
  ostream.setf(std::ios::scientific);
  ostream.precision(conf::PRECISION);
  std::size_t counter = 0;
  dipole_moment_t aveM{0.0, 0.0, 0.0};
  real_t aveM2 = 0;
  for (auto result : results.first) {
    counter += 1;

    // Get time, M, and M2.
    auto t = std::get<0>(result);
    auto M = std::get<1>(result);
    auto M2 = std::get<2>(result);
    aveM += M;
    aveM2 += M2;
    
    // Compute running averages and deviations.
    auto raveM = aveM / real_t(counter);
    auto raveM2 = aveM2 / real_t(counter);
    real_t dev = inner<real_t>(raveM, raveM) - raveM2;
    real_t rmsd = std::sqrt(dev * dev);

    // Output
    ostream << std::setw(conf::WIDTH) << t
            << std::setw(conf::WIDTH) << M
            << std::setw(conf::WIDTH) << raveM
            << std::setw(conf::WIDTH) << M2
	    << std::setw(conf::WIDTH) << raveM2
	    << std::setw(conf::WIDTH) << rmsd << std::endl;
  }
  ostream.flush();
  ostream.close();

  // Overall averages.
  aveM /= real_t(counter);
  aveM2 /= real_t(counter);

  // Probability density function for group dipole moment.
  auto fms = results.second;
  file::open_output(ostream, fnResultsFM);
  ostream.setf(std::ios::scientific);
  ostream.precision(conf::PRECISION);
  for (auto& fm: fms) {
    ostream << std::setw(conf::WIDTH) << fm.first
            << std::setw(conf::WIDTH) << fm.second * dm
            << std::setw(conf::WIDTH) << fm.first * MUUnits<real_t>::e_nm_to_D
            << std::endl;
  }
  ostream.flush();
  ostream.close();
  
  // More output.
  std::clog << std::endl;
  std::clog << "Data dipole moment M(t) written to (t is time)'"
	    << fnResultsM << "'." << std::endl;
  std::clog << "Probability P(m) of the strength m of the group dipole moment written to '"
	    << fnResultsFM << "." << std::endl;
  std::clog << "Format: t M(t) <M(t)> M2(t) "
	    << "<M2(t)> ((<M(t)>*<M(t)>-<M2(t)>)^2)^1/2" << std::endl;
  std::clog << std::endl;
  std::clog << "Average dipole moment <M>: " << aveM << std::endl;
  std::clog << "Norm of average dipole moment |<M>|: " << norm<real_t>(aveM) << std::endl;
  std::clog << "Average <M^2>: " << aveM2 << std::endl;
  real_t t = inner<real_t>(aveM, aveM) - aveM2;
  real_t rmsd = std::sqrt(t * t);
  std::clog << "( (<M>*<M> - <M2>)^2)^1/2: " << rmsd << std::endl;
  std::clog << "Dielectric constant according to the Fröhlich equation: "
	    << util::frohlich(aveM2, temperature, sm->box())
            << std::endl;

  
  
  return 0;
}
