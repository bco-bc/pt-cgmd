/*
 * File:   bond.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 15, 2019, 4:11 PM
 */

#ifndef BOND_HPP
#define BOND_HPP

#include "p-types.hpp"

namespace simploce {
    
    /**
     * A link between any two particles.
     */
    class Bond  {
    public:
    
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
        p_ptr_t getParticleOne() const;
    
        /**
         * Returns particle "2" involved in this bond.
         * @return Particle "2".
         */
        p_ptr_t getParticleTwo() const;
    
        /**
         * Is given particle part of this bond.
         * @param particle Particle.
         * @return Result.
         */
        bool contains(const p_ptr_t& particle) const;
    
    private:
    
        Bond(p_ptr_t  p1, p_ptr_t  p2) ;
    
        p_ptr_t p1_;
        p_ptr_t p2_;
    
  };
  

}


#endif /* BOND_HPP */

