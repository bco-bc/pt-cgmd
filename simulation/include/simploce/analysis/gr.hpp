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
#include "atypes.hpp"
#include "../simulation/s-properties.hpp"
#include "../simulation/bc.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/types/cvector_t.hpp"
#include "simploce/types/value_t.hpp"
#include "simploce/util/math-constants.hpp"
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <cmath>

namespace simploce {
    
    /**
     * Calculates the radial distribution function g(r) between particles of
     * given types.
     */
    template <typename P>
    class Gr : public Analyzer<P> {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = typename Analyzer<P>::p_ptr_t;
        
        /**
         * Particle group pointer type.
         */
        using pg_ptr_t = typename Analyzer<P>::pg_ptr_t;
        
        /**
         * Constructor
         * @param dr Spacing or bin size.
         * @param specName1 Particle specification name #1.
         * @param specName2 Particle specification name #2.
         * @param box Simulation box.
         * @param bc Boundary condition.
         */
        Gr(const length_t& dr,
           const std::string& specName1, 
           const std::string& specName2,
           const box_ptr_t& box,
           const bc_ptr_t& bc);
        
        void 
        perform(const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) override;
    
        /**
         * Returns results.
         * @return Radial distribution function. First of each pair is r and the second
         * is g(r).
         */
        std::vector<std::pair<real_t, real_t>> 
        results() const;
        
        /**
         * Creates an analyzer.
         * @param dr Spacing or bin size.
         * @param specName1 Particle specification name #1.
         * @param specName2 Particle specification name #2.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Analyzer.
         */
        static std::shared_ptr<Gr<P>> 
        create(const length_t& dr,
               const std::string& specName1, 
               const std::string& specName2,
               const box_ptr_t& box,
               const bc_ptr_t& bc);
    
    private:
        
        length_t dr_;
        std::string specName1_;
        std::string specName2_;
        box_ptr_t box_;
        bc_ptr_t bc_;
        
        // Helpers.
        void pp_(const std::vector<p_ptr_t>& particles);
        
        void pg_(const std::vector<p_ptr_t>& particles,
                 const std::vector<pg_ptr_t>& groups);
        
        void gg_(const std::vector<pg_ptr_t>& groups);
        
        std::size_t counter_{0};
        std::size_t nparticles1_{0};
        std::size_t nparticles2_{0}; 
        length_t rmax_{0.0};
        volume_t volume_{0.0};
        std::vector<std::size_t> hr_{};        
    };
    
    template <typename P>
    Gr<P>::Gr(const length_t& dr,
              const std::string& specName1, 
              const std::string& specName2,
              const box_ptr_t& box,
              const bc_ptr_t& bc) :
        dr_{dr}, specName1_{specName1}, specName2_{specName2}, box_{box}, bc_{bc}
    {               
        if ( specName1_.empty() || specName2_.empty() ) {
            throw std::domain_error(
                "g(r): Two particle specification names must be provided."
            );
        }
        if ( !box_ ) {
            throw std::domain_error("g(r): Missing simulation box.");
        }
        if ( !bc_ ) {
            throw std::domain_error("g(r): Missing boundary conditions.");
        }
        
        // Upper limit for g(r).
        rmax_ = util::cutoffDistance(box_);
        std::size_t nbins = rmax_() / dr_();
        hr_.resize(nbins, 0);
        volume_ = box_->volume();
    }
        
    template <typename P>
    void Gr<P>::perform(const std::vector<p_ptr_t>& all,
                        const std::vector<p_ptr_t>& free,
                        const std::vector<pg_ptr_t>& groups)
    {
        counter_ += 1;
        
        if ( counter_ == 1) {
      
            // Number of particles of the first and second specification.
            for (auto particle : all) {
                if ( particle->spec()->name() == specName1_ ) {
                    nparticles1_ += 1;
                }
                if (particle->spec()->name() == specName2_) {
                    nparticles2_ += 1;
                }
            }
            std::clog << "Number of particles of type '" << specName1_ << ": " << nparticles1_ << std::endl;
            std::clog << "Number of particles of type '" << specName2_ << ": " << nparticles2_ << std::endl;
            std::clog << "Upper limit for g(r): " << rmax_ << std::endl;
            std::clog << "Bin size: " << dr_ << std::endl;
            if ( nparticles1_ == 0 || nparticles2_ == 0 ) {
                throw std::domain_error("No such particles.");
            }
        }
        
        // For all particle pairs.
        this->pp_(all);
        
        // For all free particle pairs.
        //this->pp_(free);
        
        // For all free particle and particle groups.
        //this->pg_(free, groups);
        
        // For all groups.
        //this->gg_(groups);
        
        // Done.
    }
    
    template <typename P>
    void 
    Gr<P>::pp_(const std::vector<p_ptr_t>& particles)
    {
        static const real_t rc2 = rmax_() * rmax_();
        
        for (auto iter_i = particles.begin(); iter_i != particles.end(); ++iter_i) {
            auto pi = *iter_i;
            if ( pi->spec()->name() == specName1_ ) {
                auto ri = pi->position();
                for (auto iter_j = particles.begin(); iter_j != particles.end(); ++iter_j) {
                    if ( iter_i != iter_j) {
                        auto pj = *iter_j;
                        if ( pj->spec()->name() == specName2_ ) {
                            auto rj = pj->position();
                            auto rij = bc_->apply(ri, rj);
                            auto Rij2 = norm_square<real_t>(rij);
                            if ( Rij2 < rc2) {
                                auto Rij = std::sqrt(Rij2);
#if _DEBUG
                                util::tooClose<P>(pi, pj, std::tuple<energy_t, force_t, length_t>{});
#endif
                                std::size_t index = Rij / dr_();
                                if ( index < hr_.size() ) {
                                    hr_[index] += 1.0;
                                }
                            }
                        }
                    }
                }
            }
        }        
    }
    
    template <typename P>
    void
    Gr<P>::pg_(const std::vector<p_ptr_t>& particles,
               const std::vector<pg_ptr_t>& groups)
    {
        for (auto iter = particles.begin(); iter != particles.end(); ++iter) {
            const auto pi = *iter;
            auto ri = pi->position();
            if ( pi->spec()->name() == specName1_ ) {
                for (auto giter = groups.begin(); giter != groups.end(); ++giter) {
                    const auto g = *giter;
                    auto rg = g->position();
                    dist_vect_t r = bc_->apply(ri, rg);
                    length_t R = norm<length_t>(r);
                    if ( R() <= rmax_() ) {
                        for (auto pj : g->particles()) {
                            if ( pj->spec()->name() == specName2_ ) {
                                auto rj = pj->position();
                                dist_vect_t rij = bc_->apply(ri, rj);
                                length_t Rij = norm<length_t>(rij);
                                std::size_t index = Rij() / dr_();
                                if ( index < hr_.size() ) {
                                    hr_[index] += 1.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    template <typename P>
    void
    Gr<P>::gg_(const std::vector<pg_ptr_t>& groups)
    {
        for (auto it_i = groups.begin(); it_i != groups.end(); ++it_i) {
            auto gi = *it_i;
            auto particles_i = gi->particles();
            auto rgi = gi->position();
            for (auto it_j = groups.begin(); it_j != groups.end(); ++it_j) {
                auto gj = *it_j;
                auto rgj = gj->position();
                auto rbc = bc_->apply(rgi, rgj);
                auto Rbc = norm<length_t>(rbc);
                if ( Rbc <= rmax_() ) {
                    auto particles_j = gj->particles();
                    for (auto particle_i : particles_i) {
                        if (particle_i->spec()->name() == specName1_) {
                            auto ri = particle_i->position();
                            for (auto particle_j : particles_j) {
                                if ( particle_j->spec()->name() == specName2_ &&
                                     particle_i != particle_j) {
                                    auto rj = particle_j->position();
                                    auto rij = bc_->apply(ri, rj);
                                    auto Rij = norm<length_t>(rij);
                                    std::size_t index = Rij() / dr_();
                                    if ( index < hr_.size() ) {
                                        hr_[index] += 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    template <typename P>
    std::vector<std::pair<real_t, real_t>> Gr<P>::results() const
    {
        std::vector<std::pair<real_t, real_t>> gr{};
        
        real_t factor = 4.0 * math::constants<real_t>::PI / 3.0;

        // Number density particle of the second specification.
        real_t rho2 = nparticles2_ / volume_();

        for ( std::size_t i = 0; i != hr_.size(); ++i) {
      
            // Compute volume dV of current shell.
            real_t ri = i * dr_();
            real_t rii = ( i + 1 ) * dr_();
            real_t dV = factor * ( rii * rii * rii - ri * ri * ri );

            // Number of particles of the second specification, 
            // no correlation (ideal gas)
            real_t n = rho2 * dV;

            // Normalise, and also average over the number particles of the 
            // first specification and the number of states (observations).
            real_t g = counter_ > 0 && nparticles1_ > 0 ?
                       real_t(hr_[i]) / (n * nparticles1_ * counter_) :
                       0.0;

            // Save.
            auto pair = std::make_pair(ri, g);
            gr.push_back(pair);
        }
        
        std::clog << "Number of observations: " << counter_ << std::endl;
        
        return gr;
    }
    
    template <typename P>
    std::shared_ptr<Gr<P>> Gr<P>::create(const length_t& dr,
                                               const std::string& specName1, 
                                               const std::string& specName2,
                                               const box_ptr_t& box,
                                               const bc_ptr_t& bc)
    {
        return std::make_shared<Gr<P>>(dr, specName1, specName2, box, bc);
    }
}

#endif /* GR_HPP */

