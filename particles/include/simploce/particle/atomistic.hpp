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
 * File:   atomistic.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:19 PM
 */

#ifndef ATOMISTIC_HPP
#define ATOMISTIC_HPP

#include "particle-model.hpp"
#include "ptypes.hpp"

namespace simploce {
    
    /**
     * A physical system composed of atoms.
     */
    class Atomistic : public ParticleModel<Atom, atom_group_t> {
    public:
        
        /**
         * Constructor. No atoms, no protonation sites.
         */
        Atomistic();
        
        /**
         * Adds new atom to this physical system.
         * @param id Unique identifier.
         * @param name Name.
         * @param r Position.
         * @param spec Specification.
         * @return Newly created atom.
         */
        atom_ptr_t addAtom(std::size_t id,
                           const std::string& name,
                           const position_t& r, 
                           const atom_spec_ptr_t& spec);
        
        /**
         * Identify protonation sites.
         * @param catalog Protonation site catalog.
         */
        void protonationSites(const prot_site_catalog_ptr_t& catalog);
        
        std::size_t numberOfProtonationSites() const override;
        
        std::size_t protonationState() const override;
        
        /**
         * Returns number of atoms.
         * @return Number.
         */
        std::size_t numberOfAtoms() const;
        
        /**
         * Carries out a task with protonation sites.
         * @param task Task. Must accept protonation sites. May be a lambda expression.
         * @return Result.
         */
        template <typename R, typename T>
        R doWithProtonationSites(T task) const { return task(protonationSites_); }
        
    private:
        
        std::vector<atom_prot_site_ptr_t> protonationSites_;
    };
}

#endif /* ATOMISTIC_HPP */

