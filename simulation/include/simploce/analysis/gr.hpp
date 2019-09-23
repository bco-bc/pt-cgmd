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
#include "../simulation/bc.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/util/cvector_t.hpp"
#include "simploce/util/value_t.hpp"
#include "simploce/util/math-constants.hpp"
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

namespace simploce {
    
    /**
     * Calculates the radial distribution function g(r) for particles of
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
        
        void perform(const std::vector<p_ptr_t>& all) override;
    
        /**
         * Returns results.
         * @return Radial distribution function.
         */
        std::vector<std::pair<real_t, real_t>> results() const;
        
        /**
         * Creates an analyzer.
         * @param dr Spacing or bin size.
         * @param specName1 Particle specification name #1.
         * @param specName2 Particle specification name #2.
         * @param box Simulation box.
         * @param bc Boundary condition.
         * @return Analyzer.
         */
        static std::shared_ptr<Gr<P>> create(const length_t& dr,
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
        
        // Helpers
        std::size_t counter_{0};
        std::size_t nparticles1_{0};
        std::size_t nparticles2_{0};        
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
    }
        
    template <typename P>
    void Gr<P>::perform(const std::vector<p_ptr_t>& all)
    {
        counter_ += 1;
        
        static volume_t volume{0.0};
        static length_t rmax_{0.0};

        if ( counter_ == 1) {
            // Volume
            volume = box_->volume();
      
            // Number of particles of the first and second specification.
            for (auto particle : all) {
                if ( particle->spec()->name() == specName1_ ) {
                    nparticles1_ += 1;
                }
                if (particle->spec()->name() == specName2_) {
                    nparticles2_ += 1;
                }
            }

            // Upper limit for g(r).
            rmax_ = 0.5 * box_->size();
    
            std::size_t nbins = rmax_() / dr_();
            hr_.resize(nbins, 0);

            std::clog << "Number of particles of type '" << specName1_ << ": " << nparticles1_ << std::endl;
            std::clog << "Number of particles of type '" << specName2_ << ": " << nparticles2_ << std::endl;
            std::clog << "Upper limit for g(r): " << rmax_ << std::endl;
        }
        
        // For all particle pairs.
        for (auto iter_i = all.begin(); iter_i != all.end(); ++iter_i) {
            const P& pi = **iter_i;
            if ( pi.spec()->name() == specName1_ ) {
                position_t ri = pi.position();
                for (auto iter_j = all.begin(); iter_j != all.end(); ++iter_j) {
                    if ( iter_i != iter_j) {
                        const P& pj = **iter_j;
                        if ( pj.spec()->name() == specName2_ ) {
                            position_t rj = pj.position();
                            dist_vect_t rij = bc_->apply(ri, rj);
                            length_t Rij = norm<length_t>(rij);
                            if ( Rij() <= rmax_()) {
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
        
        // Done.
    }
    
    template <typename P>
    std::vector<std::pair<real_t, real_t>> Gr<P>::results() const
    {
        std::vector<std::pair<real_t, real_t>> gr{};
        static volume_t volume = box_->volume();
        
        real_t factor = 4.0 * MathConstants<real_t>::PI / 3.0;

        // Number density particle of the second specification.
        real_t rho2 = nparticles2_ / volume();

        for ( std::size_t i = 0; i != hr_.size(); ++i) {
      
            // Compute volume dV of current shell.
            real_t ri = i * dr_();
            real_t rii = ( i + 1 ) * dr_();
            real_t dV = factor * ( rii * rii * rii - ri * ri * ri );

            // Number of particles of the second specification, no correlation (ideal gas)
            real_t n = rho2 * dV;

            // Normalize, and also average over the number particles of the first specification as well
            // and the number of states (observations).
            real_t g = real_t(hr_[i]) / (n * nparticles1_ * counter_);

            // Save.
            auto pair = std::make_pair(ri, g);
            gr.push_back(pair);
        }
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

