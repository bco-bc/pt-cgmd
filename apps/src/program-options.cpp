/*
 * Author: André H. Juffer.
 * Created on 03/06/20222
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/util/program-options.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/param.hpp"
#include "simploce/util/file.hpp"

namespace simploce {
    namespace util {

        struct POHelper {
            std::string fnParticleSpecCatalog{"particle-spec-catalog.dat"};
            std::string fnInputParticleSystem{"in.ps"};
            std::string fnOutputParticleSystem{"out.ps"};
            std::string fnIn{"file.txt"};
            std::string fnParameters{"simulation-parameters.json"};
            std::string fnForceField{"interaction-parameters.dat"};
            spec_catalog_ptr_t catalog{};
            p_system_ptr_t particleSystem{};
            util::Logger logger{"simploce::program-options"};
        };

        static POHelper poHelper_{};

        void addStandardOptions(po::options_description &usage) {
            usage.add_options()(
                    "fn-particle-spec-catalog,s",
                    po::value<std::string>(&poHelper_.fnParticleSpecCatalog),
                    "Input file name of particle specifications. Default 'particle-spec-catalog.dat'."
            )(
                    "fn-input-particle-system,i",
                    po::value<std::string>(&poHelper_.fnInputParticleSystem),
                    "Input file name particle system. Default is 'in.ps'."
            )(
                    "fn-input-parameters,p",
                    po::value<std::string>(&poHelper_.fnParameters),
                    "Input file name parameters."
            )(
                    "fn-input-force-field,f",
                    po::value<std::string>(&poHelper_.fnForceField),
                    "Input file name interaction parameters. Default is 'interaction-parameters.dat'."
            )(
                    "fn-output-particle-system,o",
                    po::value<std::string>(&poHelper_.fnOutputParticleSystem),
                    "Output file name particle system. Default is 'out.ps'."
            )(
                    "is-coarse-grained,c",
                    "Input is a coarse-grained description. Default is an atomistic description. "
                    "No effect for input files using PDB file format."
            )(
                    "is-mesoscale,m",
                    "Simulation is at the mesoscale level, "
                    "meaning that dimensionless instead of molecular units are employed."
            )(
                    "verbose,v",
                    "Verbose"
            )(
                    "help,h",
                    "Help message"
            );
        }

        spec_catalog_ptr_t getCatalog(const po::variables_map &vm) {
            if (!poHelper_.catalog) {
                auto fn = poHelper_.fnParticleSpecCatalog;
                if (vm.count("fn-particle-spec-catalog")) {
                    fn = vm["fn-particle-spec-catalog"].as<std::string>();
                }
                poHelper_.catalog = factory::particleSpecCatalog(fn);
                poHelper_.logger.info(fn + ": Read particle specifications from this input file.");
                std::cout << *poHelper_.catalog << std::endl;
            }
            return poHelper_.catalog;
        }

        p_system_ptr_t getParticleSystem(const po::variables_map &vm, const param_ptr_t& param) {
            auto& logger = poHelper_.logger;
            if (!poHelper_.particleSystem) {
                bool isPDB{false};
                bool isCoarseGrained{false};
                if (vm.count("pdb")) {
                    isPDB = true;
                }
                if (vm.count("coarse-grained")) {
                    isCoarseGrained = true;
                    param->put("simulation.mesoscale", isCoarseGrained);
                }
                if (vm.count("is-mesoscale")) {
                    isCoarseGrained = true;
                    bool isMesoscale = true;
                    param->put<bool>("simulation.mesoscale", isMesoscale);
                }
                auto fn = poHelper_.fnInputParticleSystem;
                if (vm.count("fn-input-particle-system")) {
                    fn = vm["fn-input-particle-system"].as<std::string>();
                }
                auto catalog = util::getCatalog(vm);
                auto factory = factory::particleSystemFactory(catalog);
                if (isPDB) {
                    poHelper_.particleSystem = factory->fromPDB(fn);
                } else {
                    poHelper_.particleSystem = factory::particleSystem(fn, catalog, isCoarseGrained);
                }
                logger.info(fn + ": Read particle system from this input file.");
                auto& ps = poHelper_.particleSystem;
                logger.info(std::to_string(ps->numberOfParticles()) + ": Number of particles.");
                logger.info(std::to_string(ps->numberOfFreeParticles()) + ": Number of free particles.");
                logger.info(std::to_string(ps->numberOfParticleGroups()) + ": Number of particle groups.");
            }
            return poHelper_.particleSystem;
        }

        void writeParticleSystem(const po::variables_map& vm, const p_system_ptr_t& particleSystem) {
            auto fn = poHelper_.fnOutputParticleSystem;
            if (vm.count("fn-output-particle-system")) {
                fn = vm["fn-output-particle-system"].as<std::string>();
            }
            std::ofstream stream;
            util::open_output_file(stream, fn);
            stream << *particleSystem << std::endl;
            stream.close();
            poHelper_.logger.info("Output particle system written to '" + fn + "'.");
        }

        void verbose(const po::variables_map& vm) {
            if (vm.count("verbose")) {
                util::Logger::changeLogLevel(util::Logger::LOGTRACE);
            }
        }

        void getParameters(const po::variables_map& vm, const param_ptr_t& param) {
            std::ifstream stream;
            auto fn = poHelper_.fnParameters;
            if (vm.count("fn-input-parameters")) {
                fn = vm["fn-input-parameters"].as<std::string>();
            }
            util::open_input_file(stream, fn);
            param::read(stream, *param);
            stream.close();
        }

        ff_ptr_t getForceField(const po::variables_map& vm) {
            auto catalog = util::getCatalog(vm);
            auto fn = poHelper_.fnForceField;
            if (vm.count("fn-force-field")) {
                fn = vm["fn-force-field"].as<std::string>();
            }
            auto forceField = factory::forceField(fn, catalog);
            poHelper_.logger.info(fn + ": Read force field from this input file.");
            return std::move(forceField);
        }

        bool isMesoscale(const po::variables_map& vm) {
            return vm.count("is-mesoscale") != 0;
        }

    }
}