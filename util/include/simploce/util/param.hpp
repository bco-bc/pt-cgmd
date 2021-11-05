/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   param.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 19, 2019, 2:10 PM
 */

#ifndef PARAM_HPP
#define PARAM_HPP

#include <boost/property_tree/ptree.hpp>
#include <iostream>

namespace simploce {
    namespace param {
    
        using param_t = boost::property_tree::ptree;
    
        /**
         * Read parameters from input stream.
         * @param stream Input stream.
         * @param param Parameters,
         */
        void read(std::istream& stream, param_t& p);
    
        /**
         * Write parameters to output stream.
         * @param stream Output stream.
         * @param param Parameters.
         */
        void write(std::ostream& stream, const param_t& p);
        
        /**
         * Reads parameters from input stream.
         * @param stream Input stream.
         * @param p Parameters.
         * @return Input stream.
         */
        std::istream& operator >> (std::istream& stream, param_t& p);
        
        /**
         * Writes parameters to an output stream.
         * @param stream Output stream.
         * @param p Parameters.
         * @return Output stream.
         */
        std::ostream& operator << (std::ostream& stream, const param_t& p);
    
    }
}

#endif /* PARAM_HPP */

