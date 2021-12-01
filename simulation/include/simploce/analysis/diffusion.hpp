/*
 * Mean square deviation.
 * File:   diffusion.hpp
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on November 15, 2019, 10:45 AM
 */

#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP

#include "analyzer.hpp"
#include "a-types.hpp"
#include "simploce/simulation/s-types.hpp"
#include <memory>
#include <list>
#include <vector>
#include <utility>

namespace simploce {
    
    /**
     * Calculates group mean square deviation (msd).
     */
    template <typename P>
    class Diffusion : public analyzer<P> {
    public:
    
        /**
         * Particle pointer typeName.
         */
        using p_ptr_t = typename analyzer<P>::p_ptr_t;
        
        /**
         * Particle group pointer typeName.
         */
        using pg_ptr_t = typename analyzer<P>::pg_ptr_t;
        
        /**
         * Constructor
         * @param dt Time interval size between states.
         * @param tau Time window for msd.
         */
        Diffusion(const stime_t& dt,
                  const stime_t& tau) : dt_{dt}, tau_{tau} {}
        
        void 
        perform(const std::vector<p_ptr_t>& all,
                const std::vector<p_ptr_t>& free,
                const std::vector<pg_ptr_t>& groups) override;
        
        /**
         * Returns mean square deviation (msd), averaged over all particle groups.
         * @returns pairs of time and msd.
         */
        std::vector<std::pair<stime_t, real_t>>
        results() const;
        
        static std::shared_ptr<Diffusion<P>> create(const stime_t& dt,
                                                    const stime_t& tau);
        
    private:
        
        stime_t dt_;
        stime_t tau_;
        
        // Helpers.
        std::size_t ngroups_{};
        std::vector<real_t> msd_{};
        std::vector<std::size_t> cmsd_{};
        std::size_t counter_{0};
    };
    
    template <typename P>
    void 
    Diffusion<P>::perform(const std::vector<p_ptr_t>& all,
                          const std::vector<p_ptr_t>& free,
                          const std::vector<pg_ptr_t>& groups)
    {
        // Holds mean square deviation  over a window of length tau
        static std::size_t nmsd = 0;
        
        // Holds group position over a window of length tau.
        static std::list<std::vector<position_t>> gpos;
        
        static stime_t lastTime = 0;
        
        counter_ += 1;
        if ( counter_ == 1) {
            lastTime = counter_ * dt_;
            nmsd = tau_() / dt_();
            msd_.resize(nmsd, 0.0);
            cmsd_.resize(nmsd, 0);
            ngroups_ = groups.size();
        }
        
        // Current time
        stime_t t = counter_ * dt_;
        std::vector<position_t> positions{};
        for (auto& g : groups) {
            auto r = g->position();
            positions.push_back(r);
        }
        gpos.push_back(positions);  // At time t.
        auto tdiff = t - lastTime;
        if ( tdiff() > tau_() ) {
            gpos.pop_back();        // At time t - tau.
            lastTime += dt_;
        }
        
        std::size_t ii = 0;
        for (auto it_i = gpos.begin(); it_i != gpos.end(); ++it_i) {
            auto& r_i = *it_i; // Earlier time
            std::size_t jj = ii;
            for (auto it_j = it_i; it_j != gpos.end(); ++it_j) {
                auto& r_j = *it_j;  // Later time.
                for (std::size_t k = 0; k != r_i.size(); ++k) {
                    auto rij = r_j[k] - r_i[k];
                    auto l = jj - ii;
                    msd_[l] += norm_square<real_t>(rij);
                    cmsd_[l] += 1;
                }
                jj += 1;
            }
            ii += 1;
        }       
    }
    
    template <typename P>
    std::vector<std::pair<stime_t, real_t>>
    Diffusion<P>::results() const
    {        
        std::vector<std::pair<stime_t, real_t>> results{};
        std::size_t i = 1;
        for (auto v: msd_) {
            stime_t t = i * dt_;
            v /= (ngroups_ * cmsd_[i]);
            std::pair<stime_t, real_t> pair = std::make_pair(t, v);
            results.push_back(pair);
            i += 1;
        }        
        return results;
    }
    
    template <typename P>
    std::shared_ptr<Diffusion<P>> 
    Diffusion<P>::create(const stime_t& dt,
                         const stime_t& tau)
    {
        return std::make_shared<Diffusion<P>>(dt, tau);
    }
    
}

#endif /* DIFFUSION_HPP */

