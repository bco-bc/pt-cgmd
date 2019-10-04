#include "simploce/util/telegraph-process.hpp"
#include "simploce/util/poisson-process.hpp"
#include "simploce/util/uconf.hpp"

namespace simploce {

    using pair_t = TelegraphProcess::pair_t;

    std::vector<pair_t> 
    TelegraphProcess::generate(const stime_t& dt, 
			       const rate_t& lambda, 
			       std::size_t nValues,
			       std::size_t a)
    {
      std::vector<pair_t> result;
      int T_t = -a;
      for (std::size_t k = 0; k != nValues; ++k) {
	real_t real_k = k;
	stime_t t = real_k * dt;
	int dT = TelegraphProcess::increment(T_t, dt, lambda);
	T_t += dT;
	pair_t pair = std::make_pair(t, T_t);
	result.push_back(pair);
      }
      return result;
    }
    
    std::vector<pair_t> 
    TelegraphProcess::generate(const stime_t& dt, 
                               const rate_t& lambda_up, 
                               const rate_t& lambda_down, 
			       std::size_t nValues,
			       std::size_t a)
    {
      std::vector<pair_t> result;
      int T_t = -a;
      for (std::size_t k = 0; k != nValues; ++k) {
	real_t real_k = k;
	stime_t t = real_k * dt;
	int dT = TelegraphProcess::increment(T_t, dt, lambda_up, lambda_down);
	T_t += dT;
	pair_t pair = std::make_pair(t, T_t);
	result.push_back(pair);
      }
      return result;    
    }
    
    
    int
    TelegraphProcess::increment(int current, 
             	                const stime_t& dt, 
				const rate_t& lambda)
    {
      std::size_t dN = PoissonProcess::increment(dt, lambda);
      return -2.0 * current * dN;
    }
    
    int 
    TelegraphProcess::increment(int current, 
                                const stime_t& dt, 
				const rate_t& lambda_up, 
				const rate_t& lambda_down)
    {
      rate_t lambda = (current > 0 ? lambda_down : lambda_up);
      return TelegraphProcess::increment(current, dt, lambda);
    }
    
    std::ostream& 
    operator << (std::ostream& stream, 
                 const std::vector<pair_t>& values)
    {
      char const space = conf::SPACE;
      
      for (std::size_t k = 1; k != values.size(); ++k) {
	pair_t pair_k_1 = values[k-1];
	pair_t pair_k = values[k];
	if ( pair_k.second != pair_k_1.second ) {
	  stream << pair_k_1.first << space << pair_k.second << std::endl;
	}
	stream << pair_k.first << space << pair_k.second;
	if ( k != values.size() - 1 ) {
	  stream << std::endl;
	}
      }
      
      return stream;
    }
    
}

