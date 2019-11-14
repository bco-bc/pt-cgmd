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
 * File:   dipole-moment.hpp
 * Author: juffer
 *
 * Created on 24 October 2019, 21:21
 */

#ifndef DIPOLE_MOMENT_HPP
#define DIPOLE_MOMENT_HPP

#include "analyzer.hpp"
#include "../simulation/stypes.hpp"
#include "simploce/util/util.hpp"
#include <utility>
#include <tuple>
#include <memory>

namespace simploce {
    
    /**
     * Computes the total dipole moment M and M*M at given times, as well as
     * the probability density function f(m) of the group dipole m. 
     */
    template <typename P>
    class DipoleMoment : public Analyzer<P> {
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
         * Holds result at a single time, time, M, and M*M.
         */
        using M_result_t = std::vector<std::tuple<stime_t, dipole_moment_t, real_t>>;
        
        /**
         * Holds the probability density f(m) of the strength or size of 
         * dipole moment m of particle groups. The pair is the pair m,f(m).
         */
        using fm_result_t = std::vector<std::pair<real_t, real_t>>;
        
        /**
         * Holds the result.
         */
        using result_t = std::pair<M_result_t, fm_result_t>;
        
        /**
         * Constructor
         * @param dt Time interval between successive states.
         * @param dm Bin size of probability density f(m).
         * @param mmax Maximum value of strength m for f(m) (units are e nm). Default is
         * 10 e nm.
         * @param t0 Start time. Default is 0.
         */
        DipoleMoment(const stime_t& dt,
                     real_t dm,
                     real_t mmax = 10,
                     const stime_t& t0 = 0);
        
        void 
        perform(const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) override;
        
        /**
         * Returns results.
         * @return Results t (time), M and M*M, and m, f(m).
         */
        result_t 
        results() const;
        
        /**
         * Returns analyzer.
         * @param dt Time interval between successive states.
         * @param dm Bin size of probability density f(m).
         * @param mmax Maximum value of strength m for f(m) (units are e nm). Default is
         * 10 e nm.
         * @param t0 Start time. Default is 0.
         * @return Analyzer.
         */
        static std::shared_ptr<DipoleMoment<P>> 
        create(const stime_t& dt,
               real_t dm,
               real_t mmax = 10,
               const stime_t& t0 = 0);
        
    private:
        
        stime_t dt_;
        real_t dm_;
        real_t mmax_;
        stime_t t0_;
        M_result_t M_;
        
        // Helpers.
        std::size_t counter_{0};
        std::size_t ngroups_{0};
        std::vector<std::size_t> hm_;
    };
    
    template <typename P>
    DipoleMoment<P>::DipoleMoment(const stime_t& dt,
                                  real_t dm,
                                  real_t mmax,
                                  const stime_t& t0) : 
        dt_{dt}, dm_{dm}, mmax_{mmax}, t0_{t0}, M_{}, hm_{} 
    {                                 
        std::size_t nbins = mmax / dm_ + 1;
        hm_.resize(nbins, 0);
    }
    
    template <typename P>
    void 
    DipoleMoment<P>::perform(const std::vector<p_ptr_t>& all,
                             const std::vector<p_ptr_t>& free,
                             const std::vector<pg_ptr_t>& groups)
    {
        counter_ += 1;        
        
        // Setup if required.
        if ( counter_ == 1 ) {
            ngroups_ = groups.size();
            real_t mmax = mmax_;
            for (auto& g: groups) {
                dipole_moment_t m{0.0, 0.0, 0.0};
                for (auto& p: g->particles()) {
                    auto r = p->position();
                    auto q = p->charge();
                    for (std::size_t k = 0; k != 3; ++k) {
                        m[k] += q() * r[k];
                    }
                }
                real_t strength = norm<real_t>(m);
                if ( strength > mmax ) {
                    mmax = strength;
                }
            }
            if ( mmax > mmax_ ) {
                mmax_ = util::nint(1.25 * mmax);
                std::size_t nbins = mmax_ / dm_ + 1;
                hm_.resize(nbins, 0);    
                std::clog << "Maximum value of group dipole moment strength m " 
                          << "for probability density function f(m) set to " 
                          << mmax_ << std::endl;
            }
            std::clog << "Bin size: " << dm_ << std::endl;
        }
        
        dipole_moment_t M{0.0, 0.0, 0.0};
        for (auto p : all) {
            auto r = p->position();
            auto q = p->charge();
            for (std::size_t k = 0; k != 3; ++k) {
                M[k] += q() * r[k];
            }
        }
        real_t M2 = inner<real_t>(M, M);
        stime_t t = t0_ + counter_ * dt_;
        std::tuple<stime_t, dipole_moment_t, real_t> result = 
            std::make_tuple(t, M, M2);
        M_.push_back(result);
        
        auto last = hm_.size() - 1;
        for (auto& g: groups) {
            dipole_moment_t m{0.0, 0.0, 0.0};
            for (auto& p: g->particles()) {
                auto r = p->position();
                auto q = p->charge();
                for (std::size_t k = 0; k != 3; ++k) {
                    m[k] += q() * r[k];
                }
            }
            real_t strength = norm<real_t>(m);
            std::size_t index = strength / dm_;
            if ( index < (hm_.size() - 1) ) {
                hm_[index] += 1;
            } else {
                hm_[last] += 1;             
            }
        }
    }
    
    template <typename P>
    typename DipoleMoment<P>::result_t 
    DipoleMoment<P>::results() const
    {
        fm_result_t fm{};
        
        std::vector<real_t> hm(hm_.size(), 0.0);
        real_t tf = 0.0;
        for (std::size_t k = 0; k != hm_.size(); ++k) {
            hm[k] = hm_[k] / counter_;  // Average over states.
            tf += hm_[k];
        }
        std::size_t n = 0;
        for (auto v: hm) {
            real_t m = n * dm_;
            real_t p = v / tf;          // Probability.
            real_t f = p / dm_;         // Probability density function.
            std::pair<real_t, real_t> pair = std::make_pair(m, f);
            fm.push_back(pair);
            n += 1;
        }
        
        std::clog << "Number of observations: " << counter_ << std::endl;
        
        return std::make_pair(M_, fm);
    }
    
    template <typename P>
    std::shared_ptr<DipoleMoment<P>> 
    DipoleMoment<P>::create(const stime_t& dt,
                            real_t dm,
                            real_t mmax,
                            const stime_t& t0)
    {
        return std::make_shared<DipoleMoment<P>>(dt, dm, mmax, t0);
    }
}

#endif /* DIPOLE_MOMENT_FLUC_HPP */

