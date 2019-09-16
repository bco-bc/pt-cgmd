/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cpputil.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 12:07 PM
 */

#include "simploce/util/file.hpp"
#include <stdexcept>

namespace simploce {
    namespace file {
            
        void open_input(std::ifstream& stream, const std::string& fileName)
        {
            stream.open(fileName, std::ios_base::in);
            if ( !stream.good() ) {
                stream.close();
                throw std::runtime_error(fileName + ": Cannot open this file.");
            }
        }

        void open_output(std::ofstream& stream, const std::string& fileName)
        {
            stream.open(fileName, std::ios_base::out);
            if ( !stream.good() ) {
                stream.close();
                throw std::runtime_error(fileName + ": Cannot open this file.");
            }
        }
        
    }
}
