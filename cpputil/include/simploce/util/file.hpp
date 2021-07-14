/*
 * File:   file.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 12, 2019, 2:16 PM
 */

#ifndef FILE_HPP
#define FILE_HPP

#include <fstream>

namespace simploce {
    namespace file {
    
        /**
         * Open input file.
         * @param stream Input stream with which file is associated.
         * @param fileName Input file name.
         */
        void open_input(std::ifstream& stream, const std::string& fileName);

        /**
         * Open output file.
         * @param stream Output stream with which file is associated.
         * @param fileName Output file name.
         */
        void open_output(std::ofstream& stream, const std::string& fileName);
        
    }
}

#endif /* FILE_HPP */

