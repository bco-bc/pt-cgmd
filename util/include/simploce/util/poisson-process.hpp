#ifndef POISSON_PROCESS_HPP
#define	POISSON_PROCESS_HPP

#include "simploce/types/u-types.hpp"
#include <utility>
#include <vector>
#include <ostream>

namespace simploce {

    /**
     * Poisson process N(t), where t is time.
     * @see <a href="http://en.wikipedia.org/wiki/Poisson_process">Poisson process</a>
     */
    struct PoissonProcess {
      
      /**
       * Pair of time t and N(t).
       */
      using pair_t = std::pair<stime_t, std::size_t>;
      
      /**
       * Generates Poisson process.
       * @param dt Time step.
       * @param lambda Rate or intensity.
       * @param n Number of values (pair) to be generated..
       * @return (t,N(t)) pairs. Total number generated pairs is nValues.
       */
      static std::vector<pair_t> generate(const stime_t& dt, const rate_t& lambda, std::size_t n);
      
      /**
       * Returns increment.
       * @param dt Time interval.
       * @param lambda Rate or intensity.
       * @return Increment, either 0 or 1.
       */
      static std::size_t increment(const stime_t& dt, const rate_t& lambda);
    };
    
    /**
     * Writes Poisson process values to output stream.
     * @param stream Output stream.
     * @param values Poisson process values.
     * @return Output stream.
     */
    std::ostream& 
    operator << (std::ostream& stream, const std::vector<PoissonProcess::pair_t>& values);
    
}
  
#endif	/* POISSON_PROCESS_HPP */
  
  
