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
 * File:   pfactory.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 4:59 PM
 */

#include "simploce/particle/pfactory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-model-factory.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/file.hpp"
#include <fstream>

namespace simploce {
    namespace factory {
        
        static spec_catalog_ptr_t catalog_{};
        static particle_model_fact_ptr_t particleModelFactory_{};

        box_ptr_t cube(const length_t& side)
        {
            return std::make_shared<Cube<real_t>>(side());
        }                
        
        spec_catalog_ptr_t 
        particleSpecCatalog(const std::string& fileName)
        {
            if ( !catalog_ ) {
                std::ifstream stream;
                file::open_input(stream, fileName);
                catalog_ = factory::particleSpecCatalog(stream);
                stream.close();
            }
            return catalog_;
        }
        
        spec_catalog_ptr_t 
        particleSpecCatalog(std::istream& stream)
        {
            if ( !catalog_ ) {
                catalog_ = ParticleSpecCatalog::create(stream);
            }
            return catalog_;
        }
        
        particle_model_fact_ptr_t 
        particleModelFactory(const spec_catalog_ptr_t& catalog)
        {
            if ( !particleModelFactory_ ) {
                particleModelFactory_ = std::make_shared<ParticleModelFactory>(catalog);
            }
            return particleModelFactory_;
        }
        
        at_ptr_t atomistic()
        {
            return std::make_shared<Atomistic>();
        }
        
        cg_ptr_t coarseGrained()
        {
            return std::make_shared<CoarseGrained>();
        }
        
    }
}