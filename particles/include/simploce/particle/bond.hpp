/*
 * File:   bond.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 15, 2019, 4:11 PM
 */

#ifndef BOND_HPP
#define BOND_HPP

#include <stdexcept>
#include <memory>

namespace simploce {
    
    /**
     * A link between any two particles of any type.
     * @param P Particle type.
     */
    template <typename P>
    class Bond  {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        virtual ~Bond() {}
    
        /**
         * Creates bond between two particles.
         * @param p1 First particle, particle "1".
         * @param p2 Second particle, particle "2".
         * @return Bond.
         */
        static Bond makeBond(const p_ptr_t& p1, 
                             const p_ptr_t& p2);
    
        /**
         * Returns particle "1" involved in this bond.
         * @return Particle "1".
         */
        const p_ptr_t getParticleOne() const { return p1_; }
    
        /**
         * Returns particle "2" involved in this bond.
         * @return Particle "2".
         */
        const p_ptr_t getParticleTwo() const { return p2_; }
    
        /**
         * Is given particle part of this bond.
         * @param particle Particle.
         * @return Result.
         */
        bool contains(const p_ptr_t& particle) const;
    
    private:
    
        Bond(const p_ptr_t& p1, 
             const p_ptr_t& p2) ;
    
        p_ptr_t p1_;
        p_ptr_t p2_;
    
  };
  
  template <typename P>
  Bond<P> 
  Bond<P>::makeBond(const p_ptr_t& p1, 
                    const p_ptr_t& p2)
  {
     return Bond{p1, p2};
  }
  
  template <typename P>
  bool 
  Bond<P>::contains(const p_ptr_t& particle) const
  {
      return particle == p1_ || particle == p2_;
  }
  
  template <typename P>
  Bond<P>::Bond(const p_ptr_t& p1, 
                const p_ptr_t& p2) :
      p1_{p1}, p2_{p2}
  {      
     if ( !p1_ || !p2_ ) {
         throw std::domain_error("Bond: Two particles must be provided.");
     }
     if ( p1_ == p2_ || p1_->id() == p2_->id() ) {
         throw std::domain_error("Bond: Particles must be different.");
     }
          
  }
    
}


#endif /* BOND_HPP */

