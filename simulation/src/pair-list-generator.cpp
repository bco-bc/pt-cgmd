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
 * File:   pair-list-generator.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 11, 2019, 3:00 PM
 */

#include "simploce/simulation/pair-list-generator.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-group.hpp"

namespace simploce {
    
    static length_t RINCLUDE_{2.5};  // nm.
    
    static length_t rc2_(const box_ptr_t& box)
    {
        length_t halve = 0.5 * box->size();
        length_t rc = (RINCLUDE_() > halve() ? halve : RINCLUDE_);
        real_t rc2 = rc * rc;
        return rc2;        
    }

    template <typename P>
    struct Generator {
        
        using p_ptr_t = typename ParticlePairListGenerator<P>::p_ptr_t;
        
        using pg_ptr_t = typename ParticlePairListGenerator<P>::pg_ptr_t;
    
        using p_pair_lists_t = typename ParticlePairListGenerator<P>::p_pair_lists_t;
        
        static p_pair_lists_t create(const std::vector<p_ptr_t>& all,
                                     const std::vector<p_ptr_t>& free,
                                     const std::vector<pg_ptr_t>& groups);
    };
    
    template <typename P>
    typename Generator<P>::p_pair_lists_t
    Generator<P>::create(const std::vector<p_ptr_t>& all, 
                         const std::vector<p_ptr_t>& free, 
                         const std::vector<pg_ptr_t>& groups)
    {
        return p_pair_lists_t{};
    }
    
    ParticlePairListGenerator<Bead>::ParticlePairListGenerator(const box_ptr_t& box,
                                                               const bc_ptr_t& bc) : 
        box_{box}, bc_{bc} 
    {
    }
        
    typename ParticlePairListGenerator<Bead>::p_pair_lists_t
    ParticlePairListGenerator<Bead>::generate(const std::vector<p_ptr_t>& all,
                                              const std::vector<p_ptr_t>& free,
                                              const std::vector<pg_ptr_t>& groups) const
    {
        return Generator<Bead>::create(all, free, groups);
    }
    
}

