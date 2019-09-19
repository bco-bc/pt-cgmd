/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   param.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 19, 2019, 2:14 PM
 */

#include "simploce/util/param.hpp"
#include <boost/property_tree/json_parser.hpp>

namespace simploce {
    namespace param {
        
        void read(std::istream& stream, param_t& p)
        {
            boost::property_tree::json_parser::read_json(stream, p);
        }
        
        void write(std::ostream& stream, const param_t& p)
        {
            boost::property_tree::json_parser::write_json(stream, p);
        }
        
        std::istream& operator >> (std::istream& stream, param_t& p)
        {
            read(stream, p);
            return stream;
        }
        
        std::ostream& operator << (std::ostream& stream, const param_t& p)
        {
            write(stream, p);
            return stream;
        }
    }
}