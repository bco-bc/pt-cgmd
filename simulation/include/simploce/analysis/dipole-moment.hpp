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
 * File:   dipole-moment-fluc.hpp
 * Author: juffer
 *
 * Created on 24 October 2019, 21:21
 */

#ifndef DIPOLE_MOMENT_HPP
#define DIPOLE_MOMENT_HPP

#include "analyzer.hpp"
#include "../simulation/stypes.hpp"
#include <tuple>
#include <memory>

namespace simploce {
    
    /**
     * Computes the total dipole moment M and M*M at given times. 
     */
    template <typename P>
    class DipoleMoment : public Analyzer<P> {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = typename Analyzer<P>::p_ptr_t;
        
        /**
         * Holds result at a single time, time, M, and M*M.
         */
        using result_t = std::tuple<stime_t, dipole_moment_t, real_t>;
        
        /**
         * Holds results for all times.
         */
        using results_t = std::vector<result_t>;
        
        /**
         * Constructor
         * @param dt Time interval between successive states.
         * @param t0 Start time. Default is 0.
         */
        DipoleMoment(const stime_t& dt,
                     const stime_t& t0 = 0);
        
        void perform(const std::vector<p_ptr_t>& all) override;
        
        /**
         * Returns results.
         * @return Results t (time), M and M*M.
         */
        results_t results() const { return results_; }
        
        /**
         * Returns analyzer.
         * @param dt Time interval between successive states.
         * @param t0 Start time. Default is 0.
         * @return Analyzer.
         */
        static std::shared_ptr<DipoleMoment<P>> create(const stime_t& dt,
                                                       const stime_t& t0 = 0);
        
    private:
        
        stime_t dt_;
        stime_t t0_;
        results_t results_;
        
        // Helpers.
        std::size_t counter_{0};
        
        
    };
    
    template <typename P>
    DipoleMoment<P>::DipoleMoment(const stime_t& dt,
                                  const stime_t& t0) : dt_{dt}, t0_{t0}, results_{}
    {                                  
    }
    
    template <typename P>
    void DipoleMoment<P>::perform(const std::vector<p_ptr_t>& all)
    {
        counter_ += 1;
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
        result_t result = std::make_tuple(t, M, M2);
        results_.push_back(result);        
    }
    
    template <typename P>
    std::shared_ptr<DipoleMoment<P>> DipoleMoment<P>::create(const stime_t& dt,
                                                             const stime_t& t0)
    {
        return std::make_shared<DipoleMoment<P>>(dt, t0);
    }
}

#endif /* DIPOLE_MOMENT_FLUC_HPP */

