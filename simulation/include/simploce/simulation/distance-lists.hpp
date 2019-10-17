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
 * File:   distance-lists.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 17, 2019, 12:13 PM
 */

#ifndef DISTANCE_LISTS_HPP
#define DISTANCE_LISTS_HPP

#include "pair-list-generator.hpp"
#include "simploce/particle/bead.hpp"
#include "stypes.hpp"

namespace simploce {
    
    /**
     * Creates particle pair lists based on distances between particles.
     * @param P Particle type.
     */
    template <typename P>
    class DistanceLists;
    
    /**
     * Specialization for atoms
     */
    template <>
    class DistanceLists<Atom> : public ParticlePairListGenerator<Atom> {
    public:
        
        using p_ptr_t = typename ParticlePairListGenerator<Atom>::p_ptr_t;
        using pg_ptr_t = typename ParticlePairListGenerator<Atom>::pg_ptr_t;
        using p_pair_t = typename ParticlePairListGenerator<Atom>::p_pair_t;
        using p_pair_list_t = typename ParticlePairListGenerator<Atom>::p_pair_list_t;
                
        DistanceLists(const box_ptr_t& box,
                      const bc_ptr_t& bc);
        
        std::vector<p_pair_list_t> 
        generate(const std::vector<p_ptr_t>& all,
                 const std::vector<p_ptr_t>& free,
                 const std::vector<pg_ptr_t>& groups) const override;
        
    private:
        
        box_ptr_t box_;
        bc_ptr_t bc_;
                
    };
    
    
    /**
     * Specialization for beads.
     */
    template <>
    class DistanceLists<Bead> : public ParticlePairListGenerator<Bead> {
    public:
        
        using p_ptr_t = typename ParticlePairListGenerator<Bead>::p_ptr_t;
        using pg_ptr_t = typename ParticlePairListGenerator<Bead>::pg_ptr_t;
        using p_pair_t = typename ParticlePairListGenerator<Bead>::p_pair_t;
        using p_pair_list_t = typename ParticlePairListGenerator<Bead>::p_pair_list_t;
                
        DistanceLists(const box_ptr_t& box,
                      const bc_ptr_t& bc);
        
        std::vector<p_pair_list_t> 
        generate(const std::vector<p_ptr_t>& all,
                 const std::vector<p_ptr_t>& free,
                 const std::vector<pg_ptr_t>& groups) const override;
        
    private:
        
        box_ptr_t box_;
        bc_ptr_t bc_;
                
    };
}


#endif /* DISTANCE_LISTS_HPP */

