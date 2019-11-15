/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   utypes.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 4, 2019, 3:38 PM
 */

#ifndef UTYPES_HPP
#define UTYPES_HPP

#include "value_t.hpp"
#include "cvector_t.hpp"
#include "cube.hpp"

namespace simploce {
    
    /**
     * Real numbers type.
     */
    using real_t = double;
    
    /**
     * Time type.
     */
    using stime_t = value_t<real_t, -1>;
    
    /**
     * Rate type.
     */
    using rate_t = value_t<real_t, -2>;
    
    
}

#endif /* UTYPES_HPP */

