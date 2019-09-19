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
 * File:   atom.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 2:58 PM
 */

#include "simploce/particle/atom.hpp"
#include <boost/lexical_cast.hpp>
#include <stdexcept>

namespace simploce {
    
    Atom::Atom(std::size_t index, const std::string& name, const atom_spec_ptr_t& spec) :
        Particle(index, name, spec)
    {        
    }
    
    Atom::~Atom()
    {        
    }
    
    atom_ptr_t Atom::create(std::size_t index, 
                            const std::string &name, 
                            const atom_spec_ptr_t &spec)
    {
        return atom_ptr_t(new Atom(index, name, spec));
    }
}