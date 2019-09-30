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
 * File:   newfile.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 4:59 PM
 */

#include "simploce/particle/pfactory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/util/file.hpp"
#include <fstream>

namespace simploce {
    namespace factory {
        
        spec_catalog_ptr_t 
        particleSpecCatalog(const std::string& fileName)
        {
            std::ifstream stream;
            file::open_input(stream, fileName);
            spec_catalog_ptr_t catalog = factory::particleSpecCatalog(stream);
            stream.close();
            return catalog;
        }
        
        spec_catalog_ptr_t 
        particleSpecCatalog(std::istream& stream)
        {
            return ParticleSpecCatalog::create(stream);
        }
    }
}