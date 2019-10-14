/*
 * The MIT License
 *
 * Copyright 2019 juffer.
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
 * File:   cell-lists.hpp
 * Author: juffer
 *
 * Created on 11 October 2019, 16:29
 */

#ifndef CELL_LISTS_HPP
#define CELL_LISTS_HPP

#include "pair-list-generator.hpp"
#include "stypes.hpp"
#include "simploce/particle/bead.hpp"

namespace simploce {
    
    /**
     * Finds all particles pairs within a given cut-off distance of each other 
     * in molecular dynamics simulations.
     * @see <a href="https://en.wikipedia.org/wiki/Cell_lists">Cell lists at Wikipedia</a>
     */
    template <typename P>
    class CellLists;
    
    /**
     * Specialisation for beads.
     */
    template <>
    class CellLists<Bead> : public ParticlePairListGenerator<Bead> {
    public:
        
        using p_ptr_t = typename ParticlePairListGenerator<Bead>::p_ptr_t;
        using pg_ptr_t = typename ParticlePairListGenerator<Bead>::pg_ptr_t;
        using p_pair_list_t = typename ParticlePairListGenerator<Bead>::p_pair_list_t;
        
        CellLists(const box_ptr_t& box,
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

#endif /* CELL_LISTS_HPP */

