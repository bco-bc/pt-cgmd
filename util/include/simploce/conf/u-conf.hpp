/*
 * File:   conf.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 3, 2019, 11:28 AM
 */

#ifndef U_CONF_HPP
#define U_CONF_HPP

namespace simploce {
    /**
     * Some universal settings.
     */
    namespace conf {
        
        /**
         * Precision of numbers in output streams.
         */
        const int PRECISION = 8;

        /**
         * Width of real number in output streams.
         */
        const int REAL_WIDTH = 16;
        
        /**
         * Width of integer numbers in output streams.
         */
        const int INTEGER_WIDTH = 10;

        /**
         * Width of identifiers in output streams.
         */
        const int ID_WIDTH = 17;
        
        /**
         * Space symbol in output streams.
         */
        const char SPACE = ' ';


    }
}

#endif /* U_CONF_HPP */

