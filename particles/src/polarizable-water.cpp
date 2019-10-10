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
#include "simploce/particle/bead.hpp"
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
        std::vector<p_ptr_t>& particles = this->all();

        std::random_device rd;  
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, groups.size()-1);
 
        // Select a group.
        std::size_t k = dis(gen);
        auto& group = groups[k];
        
        // Erase the group, keep group particles.
        position_t r = group->position();
        auto groupParticles = group->particles();
        groups.erase(groups.begin() + k);
        
        // Remove group particles.
        std::vector<std::vector<p_ptr_t>::iterator> iters;
        for (auto iter = groupParticles.begin(); iter != groupParticles.end(); ++iter) {
            auto bead = *iter;
            auto it = particles.begin() + bead->index();
            iters.push_back(it);
        }
        for (auto iter = iters.begin(); iter != iters.end(); ++iter) {
            particles.erase(*iter);
        }
        
        return r;
    }
    
}