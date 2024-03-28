/*
 * Author: André H. Juffer.
 * Created on 19/11/2021, 14:06.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef APPS_PROGRAM_OPTIONS_HPP
#define APPS_PROGRAM_OPTIONS_HPP

#include "simploce/types/s-types.hpp"
#include "simploce/particle/p-factory.hpp"
#include <boost/program_options.hpp>
#include <string>

namespace po = boost::program_options;

namespace simploce {
    namespace util {

        /**
         * Adds the standard and/or commonly options for an application.
         * @param usage Usage description.
         */
        void addStandardOptions(po::options_description& usage);

        /**
         * Returns particle specifications catalog, read from an input file.
         * @param vm Programs options and values.
         * @param fileName Default file name, if not available from program options and values,
         * @return Particle specifications catalog.
         */
        spec_catalog_ptr_t getCatalog(const po::variables_map& vm);

        /**
         * Returns particle system, read from an input file.
         * @param vm Program options and values.
         * @param Parameters.
         * @return Particle system.
         */
        p_system_ptr_t getParticleSystem(const po::variables_map& vm, const param_ptr_t& param);

        /**
         * Writes particle system.
         * @param vm Program options and values.
         * @param particleSystem Particle system.
         */
        void writeParticleSystem(const po::variables_map& vm, const p_system_ptr_t& particleSystem);

        /**
         * Sets the debugging logging level. This may result in a large output log file, if a verbose log
         * level is requested.
         * @param vm Program options and values.
         */
        void verbose(const po::variables_map& vm);

        /**
         * Returns parameters.
         * @param vm Program options and values.
         * @param param Existing parameters, may not contain anything, but must not be a nullptr. On return, its
         * elements may be overwritten.
         */
        void getParameters(const po::variables_map& vm, const param_ptr_t& param);

        /**
         * Returns force field, read from an input file.
         * @param vm Program options and values.
         * @return Force field.
         */
        ff_ptr_t getForceField(const po::variables_map& vm);

        /**
         * Mesoscale?
         * @param vm Program options and values.
         * @return Result.
         */
        bool isMesoscale(const po::variables_map& vm);

    }
}

#endif //APPS_PROGRAM_OPTIONS_HPP
