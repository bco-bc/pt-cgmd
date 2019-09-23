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
 * File:   gr.hpp
 * Author: juffer
 *
 * Created on 22 September 2019, 15:22
 */

#ifndef GR_HPP
#define GR_HPP

#include "analyzer.hpp"
#include <memory>
#include <string>

namespace simploce {
    
    /**
     * Calculates the radial distribution function g(r) for particles of
     * given types.
     */
    template <typename P>
    class Gr : public Analyzer<P> {
    public:
        
        /**
         * Constructor
         * @param specName1 Particle specification name #1.
         * @param specName2 Particle specification name #2.
         */
        Gr(const std::string& specName1, 
           const std::string& specName2);
        
        void perform(const std::vector<p_ptr_t>& all,
                     const box_ptr_t& box,
                     const bc_ptr_t& bc) const override;
    
        /**
         * Returns results.
         * @return Radial distribution function.
         */
        std::vector<std::pair<real_t, real_t>> results() const;
    
    private:
        
        std::string specName1_;
        std::string specName2_;
        
        std::size_t counter_;
        std::vector<std::pair<real_t, real_t>> gr_;        
    };
    
    template <typename P>
    Gr::Gr(const std::string& specName1, 
           const std::string& specName2) :
        specName1_{specName1}, specName2_{specName2}, counter_{0}, gr_{}
    {               
    }
        
    template <typename P>
    void Gr::perform(const std::vector<p_ptr_t> &particles,
                     const box_ptr_t& box,
                     const bc_ptr_t& bc)
    {
        
    }
}

#endif /* GR_HPP */

