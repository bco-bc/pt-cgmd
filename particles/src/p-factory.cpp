/*
 * File:   p-factory.cpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on September 6, 2019, 4:59 PM
 */

#include "simploce/particle/p-factory.hpp"
#include "simploce/particle/particle-spec-catalog.hpp"
#include "simploce/particle/particle-system-factory.hpp"
#include "simploce/particle/atomistic.hpp"
#include "simploce/particle/coarse-grained.hpp"
#include "simploce/util/file.hpp"
#include <fstream>

namespace simploce {
    namespace factory {
        
        static spec_catalog_ptr_t catalog_{};
        static p_system_fact_ptr_t particleModelFactory_{};

        box_ptr_t box(const length_t& side)
        {
            return std::make_shared<Cube<real_t>>(side());
        }                
        
        spec_catalog_ptr_t 
        particleSpecCatalog(const std::string& fileName)
        {
            if ( !catalog_ ) {
                std::ifstream stream;
                util::open_input_file(stream, fileName);
                catalog_ = factory::particleSpecCatalog(stream);
                stream.close();
            }
            return catalog_;
        }
        
        spec_catalog_ptr_t 
        particleSpecCatalog(std::istream& stream)
        {
            if ( !catalog_ ) {
                catalog_ = ParticleSpecCatalog::obtainFrom(stream);
            }
            return catalog_;
        }
        
        p_system_fact_ptr_t
        particleModelFactory(const spec_catalog_ptr_t& catalog)
        {
            if ( !particleModelFactory_ ) {
                particleModelFactory_ = std::make_shared<ParticleSystemFactory>(catalog);
            }
            return particleModelFactory_;
        }
        
        at_sys_ptr_t atomistic()
        {
            auto atomistic = std::make_shared<Atomistic>();
            auto box = factory::box(0.0);
            atomistic->box(box);
            return atomistic;
        }
        
        cg_sys_ptr_t coarseGrained()
        {
            auto cg = std::make_shared<CoarseGrained>();
            auto box = factory::box(0.0);
            cg->box(box);
            return cg;
        }
        
    }
}