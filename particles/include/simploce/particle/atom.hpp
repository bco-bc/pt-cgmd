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
 * File:   atom.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on July 30, 2019, 4:46 PM
 */

#ifndef ATOM_HPP
#define ATOM_HPP

#include "particle.hpp"
#include "ptypes.hpp"
#include <string>

namespace simploce {
    
    /**
     * The smallest constituent unit of ordinary matter that has the 
     * properties of a chemical element.
     */
    class Atom: public Particle {
    public:
        
        // Noncopyable.
        Atom(const Atom&) = delete;
        Atom& operator = (const Atom&) = delete;
        
        ~Atom();
                
    private:
        
        friend class Atomistic;
        
        static atom_ptr_t create(std::size_t index, 
                                 const std::string &name, 
                                 const atom_spec_ptr_t &spec);
        
        Atom(std::size_t index, 
             const std::string &name,
             const atom_spec_ptr_t &spec);        
        
    };
}

#endif /* ATOM_HPP */

