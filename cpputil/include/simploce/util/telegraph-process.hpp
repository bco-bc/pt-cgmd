#ifndef TELEGRAPH_PROCESS_HPP
#define TELEGRAPH_PROCESS_HPP

#include "utypes.hpp"
#include <utility>
#include <vector>
#include <ostream>

namespace simploce {
        
    /**
     * Telegraph process T(t), where t is time. Increments are generated from the 
     * increment of a Poisson process. The process jumps between two integral 
     * values +a (up, right) and -a (down, left). Here, a > 0 and is a natural number.
     * @see <a href="http://en.wikipedia.org/wiki/Telegraph_process">Telegraph process</a>
     */
    struct TelegraphProcess {
      
      /**
       * A pair of t and T(t).
       */
      using pair_t = std::pair<stime_t, int>;
      
      /**
       * Generates symmetric Telegraph process jumping between -a and +a, a > 0.
       * @param dt Time step.
       * @param lambda Rate.
       * @param nValues Requested number of values.
       * @param a Value of a. Should be a > 0.
       * @return Tuple of (t,T(t)) pairs. Initial value T(0) = -a.
       */
      static std::vector<pair_t> generate(const stime_t& dt, 
					  const rate_t& lambda, 
					  std::size_t nValues,
					  std::size_t a);
      
      /**
       * Generates asymmetric Telegraph process jumping between -a and +a, a > 0.
       * @param dt Time step.
       * @param lambda_up Rate for going from -a to +a.
       * @param lambda_down Rate for going from +a to -a.
       * @param nValues Requested number of values.
       * @param a Value of a. Should be a > 0.
       * @return Tuple of (t,T(t)) pairs.
       */
      static std::vector<pair_t> generate(const stime_t& dt, 
					  const rate_t& lambda_up, 
					  const rate_t& lambda_down, 
					  std::size_t nValues,
					  std::size_t a);
      
      /**
       * Returns increment for a symmetric telegraph process for jumping between 
       * -a and +a,
       * with a = | current |.
       * @param current Current state value
       * @param dt Time interval.
       * @param lambda Jump rate.
       * @return Increment, either 0, +2 * | current |, or -2 * | current |.
       */
      static int increment(int current, 
		           const stime_t& dt, 
			   const rate_t& lambda);
      
      /**
       * Returns increment for a asymmetric telegraph for jumping between -a and +a, 
       * with a = | current |.
       * @param current Current state value.
       * @param dt Time interval.
       * @param rate_up Rate for jumping from -a to +a.
       * @param rate_down Rate for jumping from +a to -a.
       * @return Increment, either 0, +2 * |current|, or -2 * |current|.
       */
      static int increment(int current, 
		           const stime_t& dt, 
			   const rate_t& rate_up, 
			   const rate_t& rate_down);
    };
    
    /**
     * Writes Telegraph process values to output stream.
     * @param stream Output stream.
     * @param values Telegraph process values.
     * @return Output stream.
     */
    std::ostream& 
    operator << (std::ostream& stream, 
                 const std::vector<TelegraphProcess::pair_t>& values);
    
}

#endif	/* TELEGRAPH_PROCESS_HPP */

