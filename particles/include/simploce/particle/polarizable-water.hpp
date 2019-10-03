/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   polarizable-water.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 2:21 PM
 */

#ifndef POLARIZABLE_WATER_HPP
#define POLARIZABLE_WATER_HPP

#include "coarse-grained.hpp"

namespace simploce {
    
    /**
     * Coarse grained particle model for polarizable water.
     */
    class PolarizableWater : public CoarseGrained {
    public:
        
        PolarizableWater();
        
    private:
        
        friend class ParticleModelFactory;
        
        position_t removeGroup_() override;
        
    };
}

#endif /* POLARIZABLE_WATER_HPP */

