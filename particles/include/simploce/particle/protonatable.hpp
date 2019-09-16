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
 * File:   protonatable.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 11:37 AM
 */

#ifndef PROTONATABLE_HPP
#define PROTONATABLE_HPP

#include <cstdlib>

namespace simploce {
    
    /**
     * Capable of binding and releasing protons.
     */
    struct Protonatable {
        
        ~Protonatable() {}
        
        /**
         * Binds one more proton.
         */
        virtual void protonate() = 0;
        
        /**
         * Releases one proton.
         */
        virtual void deprotonate() = 0;
        
        /**
         * Inquires whether this entity is already protonated.        
         * @return Result.
         */
        virtual bool isProtonated() const = 0;
        
        /**
         * Returns number of bound protons.
         * @return Number, always >= 0.
         */
        virtual std::size_t protonationState() const = 0;
    };
}

#endif /* PROTONATABLE_HPP */

