/*
 * The MIT License
 *
 * Copyright 2019 André H. Juffer, Biocenter Oulu
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
 * File:   protonation-site.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 5, 2019, 4:43 PM
 */

#ifndef PROTONATION_SITE_HPP
#define PROTONATION_SITE_HPP

#include "protonatable.hpp"
#include "particle.hpp"
#include <stdexcept>
#include <vector>
#include <memory>
#include <cassert>

namespace simploce {
    
    /**
     * A protonatable entity composed of two or more particles. (De)protonation 
     * results in a redistribution of charge over the constituting particles. 
     * A maximum of one proton can be bound or released.
     * @param T Particle type.
     * @param S Particle specification type.
     */
    template <typename P, typename S>
    class ProtonationSite: public Protonatable {
    public:
        
        /**
         * Particle pointer type.
         */
        using p_ptr_t = std::shared_ptr<P>;
        
        /**
         * Particle specification type.
         */
        using ps_ptr_t = std::shared_ptr<S>;        
        
        /**
         * Constructor. Creates a protonation site, which initially will be deprotonated.
         * @param name Name, e.g. COOH.
         * @param particles Constituting particles.
         * @param deprotonated Specification for the deprotonated state. The order of 
         * specifications must be according to the constituting particles. Not checked!
         * @param protonated Specifications for the protonated state. The order of 
         * specifications must be according to the constituting particles. Not checked!
         */
        ProtonationSite(const std::string& name,
                        const std::vector<p_ptr_t> &particles,
                        const std::vector<ps_ptr_t> &deprotonated,
                        const std::vector<ps_ptr_t> &protonated);
        
        /**
         * Returns site name.
         * @return Name.
         */
        std::string name() const;
    
        void protonate() override;
        
        void deprotonate() override;
        
        bool isProtonated() const override;
        
        std::size_t protonationState() const override;
        
        /**
         * Returns total charge value.
         * @return Value.
         */
        charge_t charge() const;
        
    private:
        
        std::string name_;
        std::vector<p_ptr_t> particles_;
        std::vector<ps_ptr_t> deprotonated_;
        std::vector<ps_ptr_t> protonated_;
        bool proton_;
    };
    
    template <typename T, typename S>
    ProtonationSite<T,S>::ProtonationSite(const std::string& name,
                                          const std::vector<p_ptr_t>& particles,
                                          const std::vector<ps_ptr_t>& deprotonated,
                                          const std::vector<ps_ptr_t>& protonated) :
        name_(name), particles_{particles}, deprotonated_{deprotonated}, 
        protonated_{protonated}, proton_{true}
    {
        if ( name_.empty() ) {
            throw std::domain_error(
                "ProtonationSite: Name must be provided."
            );
        }
        std::size_t size = particles_.size();
        if ( size <= 1 ) {
            throw std::domain_error(
                "ProtonationSite: Must consist of two or more particles."
            );
        }
        if ( size != deprotonated_.size() || size != protonated_.size() ) {
            throw std::domain_error(
                "ProtonationSite: Number of particles is not equal to number of "
                "particle specifications."
            );
        }   
        
        // Initial state is the deprotonated.
        this->deprotonate();
    }
        
    template <typename T, typename S>
    std::string ProtonationSite<T,S>::name() const 
    {
        return name_;
    }
    
    template <typename T, typename S>
    void ProtonationSite<T,S>::protonate()
    {
        assert( !this->isProtonated() );
        for (std::size_t index = 0; index != particles_.size(); ++index) {
            auto p = particles_[index];
            auto spec = protonated_[index];
            p->reset_(spec);
        }
        proton_ = true;
    }
    
    template <typename T, typename S>
    void ProtonationSite<T,S>::deprotonate()
    {
        assert( this->isProtonated() );
        for (std::size_t index = 0; index != particles_.size(); ++index) {
            auto p = particles_[index];
            auto spec = deprotonated_[index];
            p->reset_(spec);
        }
        proton_ = false;
    }
    
    template <typename T, typename S>
    bool ProtonationSite<T,S>::isProtonated() const
    {
        return proton_;
    }
    
    template <typename T, typename S>
    std::size_t ProtonationSite<T,S>::protonationState() const
    {
        return (this->isProtonated() ? 1 : 0);
    }
    
    template <typename T, typename S>
    charge_t ProtonationSite<T,S>::charge() const
    {
        const std::vector<ps_ptr_t>& state = 
            this->isProtonated() ? protonated_ : deprotonated_;
        charge_t charge{0};
        for (auto spec : state) {
            charge += spec->charge();
        }
        return charge;
    }
     
}


#endif /* PROTONATION_SITE_HPP */

