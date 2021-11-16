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
#include "s-types.hpp"
#include "simploce/particle/bead.hpp"
#include "simploce/particle/atom.hpp"

namespace simploce {
    
    /**
     * Creates particle pair lists based on location of particles in cells.
     * @see <a href="https://en.wikipedia.org/wiki/Cell_lists">Cell lists at Wikipedia</a>
     */
    template <typename P>
    class CellLists;
    
    /**
     * Specialization for atoms.
     */
    template <>
    class CellLists<Atom> : public pair_lists_generator<Atom> {
    public:
        
        
        CellLists(const box_ptr_t& box,
                  const bc_ptr_t& bc);
        
        PairLists<Atom>
        generate(const std::vector<atom_ptr_t>& all,
                 const std::vector<atom_ptr_t>& free,
                 const std::vector<atom_group_ptr_t>& groups) const override;
        
        /**
         * Sets and override cell side length. Default is 0.5 of cutoff distance.
         * @param sideLength Side length.
         */
        void sideLength(const length_t& sideLength);
        
    private:
                        
        box_ptr_t box_;
        bc_ptr_t bc_;
        length_t sideLength_;
    };
    
    
    /**
     * Specialization for beads.
     */
    template <>
    class CellLists<Bead> : public pair_lists_generator<Bead> {
    public:
        
        CellLists(const box_ptr_t& box,
                  const bc_ptr_t& bc);
        
        PairLists<Bead> 
        generate(const std::vector<bead_ptr_t>& all,
                 const std::vector<bead_ptr_t>& free,
                 const std::vector<bead_group_ptr_t>& groups) const override;
        /**
         * Sets and override cell side length. Default is 0.5 of cutoff distance.
         * @param sideLength Side length.
         */
        void sideLength(const length_t& sideLength);
        
    private:
                        
        box_ptr_t box_;
        bc_ptr_t bc_;
        length_t sideLength_;
    };
}

#endif /* CELL_LISTS_HPP */

