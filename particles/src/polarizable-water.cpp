/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   polarizable-water.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on October 3, 2019, 3:01 PM
 */

#include "simploce/particle/polarizable-water.hpp"
#include "simploce/particle/particle-spec.hpp"
#include <random>

namespace simploce {
    
    PolarizableWater::PolarizableWater() : 
        CoarseGrained{}
    {        
    }
    
    position_t 
    PolarizableWater::removeGroup_()
    {
        std::vector<pg_ptr_t>& groups = this->groups();

        std::random_device rd;  
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, groups.size()-1);
 
        std::size_t k = dis(gen);
        position_t r = groups[k]->position();
        auto iter = groups.begin() + k;
        groups.erase(iter);
        
        return r;
    }
    
}