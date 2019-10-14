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
#include "simploce/simulation/cell-lists.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/atom.hpp"
#include "simploce/particle/particle-group.hpp"

namespace simploce {
    
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
        static CellLists<Bead> cellLists{box_, bc_};
        
        return cellLists.(all, free, groups);
    }
    
}

