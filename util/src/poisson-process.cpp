/*
 * File:   poisson-process.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 19, 2019, 2:14 PM
 */

#include "simploce/util/poisson-process.hpp"
#include "simploce/conf/u-conf.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <numeric>
#include <string>
#include <random>
#include <iostream>

namespace simploce {

    using pair_t = PoissonProcess::pair_t;

    std::vector<pair_t> 
    PoissonProcess::generate(const stime_t&  dt, const rate_t& lambda, std::size_t n) {
        std::vector<pair_t> result;
        std::size_t N_t = 0;
        for (std::size_t k = 0; k != n; ++k ) {
	        real_t real_k = k;
	        stime_t t = real_k * dt;
	        std::size_t dN = PoissonProcess::increment(dt, lambda);
	        N_t += dN;
	        pair_t pair = std::make_pair(t, N_t);
	        result.push_back(pair);
      }
      return result;
    }
    
    std::size_t 
    PoissonProcess::increment(const stime_t& dt, const rate_t& lambda) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<real_t> dis(0.0, 1.0);
      
        real_t v = dis(gen);
        real_t exp = std::exp( -lambda() * dt() );
        std::size_t result =  ( v > exp ? 1 : 0 );
        return result;
    }
    
    std::ostream& 
    operator << (std::ostream& stream, const std::vector<pair_t>& values) {
        for ( std::size_t k = 1; k != values.size(); ++k) {
	        pair_t pair_k_1 = values[k-1];
	        pair_t pair_k = values[k];
	        if ( pair_k.second > pair_k_1.second ) {
	            stream << std::setw(conf::REAL_WIDTH) << pair_k_1.first;
                stream << conf::SPACE << std::setw(conf::REAL_WIDTH) << pair_k.second << std::endl;
	        }
	    stream << pair_k.first << conf::SPACE << std::setw(conf::REAL_WIDTH) << pair_k.second;
	    if ( k != values.size() - 1 ) {
	        stream << std::endl;
	    }
      }
      return stream;
    }

}
