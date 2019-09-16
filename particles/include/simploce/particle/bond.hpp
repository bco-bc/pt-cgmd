/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
        using particle_ptr_t = std::shared_ptr<P>;
        
        virtual ~Bond() {}
    
        /**
         * Creates bond between two particles.
         * @param p1 First particle, particle "1".
         * @param p2 Second particle, particle "2".
         * @return Bond.
         */
        static Bond makeBond(const particle_ptr_t& p1, const particle_ptr_t& p2);
    
        /**
         * Returns particle "1" involved in this bond.
         * @return Particle "1".
         */
        const particle_ptr_t getParticleOne() const { return p1_; }
    
        /**
         * Returns particle "2" involved in this bond.
         * @return Particle "2".
         */
        const particle_ptr_t getParticleTwo() const { return p2_; }
    
        /**
         * Is given particle part of this bond.
         * @param particle Particle.
         * @return Result.
         */
        bool contains(const particle_ptr_t& particle) const;
    
    private:
    
        Bond(const particle_ptr_t& p1, const particle_ptr_t& p2) ;
    
        particle_ptr_t p1_;
        particle_ptr_t p2_;
    
  };
  
  template <typename P>
  Bond<P> Bond<P>::makeBond(const particle_ptr_t& p1, const particle_ptr_t& p2)
  {
     if ( !p1 || !p2 ) {
         throw std::domain_error("Bond: Two particles must be provided.");
     }
     return Bond{p1, p2};
  }
  
  template <typename P>
  bool Bond<P>::contains(const particle_ptr_t& particle) const
  {
      return particle == p1_ || particle == p2_;
  }
  
  template <typename P>
  Bond<P>::Bond(const particle_ptr_t& p1, const particle_ptr_t& p2) :
      p1_{p1}, p2_{p2}
  {      
  }
    
}


#endif /* BOND_HPP */

