/*
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

        void write(std::ostream& stream, const param_ptr_t& p) {
            write(stream, *p);
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