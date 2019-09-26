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
 * File:   coarse-grained.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 7, 2019, 2:34 PM
 */

#ifndef COARSE_GRAINED_HPP
#define COARSE_GRAINED_HPP

#include "particle-model.hpp"
#include "ptypes.hpp"
#include <utility>

namespace simploce {
    
    /**
     * A physical system composed of regular and protonatable beads.
     */
    class CoarseGrained : public ParticleModel<Bead, bead_group_t> {
    public:
        
        /**
         * Constructor. Empty particle model.
         */
        CoarseGrained();
        
        /**
         * Adds a new bead to this physical system.
         * @param id Unique identifier.
         * @param name Name.
         * @param r Position.
         * @param spec Specification.
         */
        bead_ptr_t addBead(std::size_t id,
                           const std::string& name, 
                           const position_t& r, 
                           const bead_spec_ptr_t& spec);
        
        /**
         * Adds new discrete protonatable bead to this coarse grained particle model.
         * @param id Unique identifier.
         * @param name Name.
         * @param r Position.
         * @param protonationState Number of bound protons. Must be >= 0.
         * @param spec Specification.
         */
        dprot_bead_ptr_t addProtonatableBead(std::size_t id,
                                             const std::string& name, 
                                             const position_t& r,
                                             std::size_t protonationState,
                                             const bead_spec_ptr_t& spec);
        
        /**
         * Adds new protonating bead to this coarse grained particle model.
         * @param id Unique identifier.
         * @param name Name.
         * @param r Position.
         * @param protonationState Number of bound protons. Must be >= 0.
         * @param spec Specification.
         */
        cprot_bead_ptr_t addProtonatingBead(std::size_t id,
                                            const std::string& name, 
                                            const position_t& r,
                                            std::size_t protonationState,
                                            const bead_spec_ptr_t& spec);
                        
        /**
         * Adds a bead group with bonds to this physical system.
         * @param beads Beads. Must be already in this particle model.
         * @param bbonds Bonds, given as pairs of bead identifiers.
         * @return Added bead group.
         */
        bead_group_ptr_t addBeadGroup(const std::vector<bead_ptr_t>& beads, 
                                      const std::vector<id_pair_t>& bbonds);
        
        /**
         * Returns the total number of discrete and continuous protonatables.
         * @return 
         */
        std::size_t numberOfProtonationSites() const override;
        
        std::size_t protonationState() const override;
        
        /**
         * Reads a coarse grained particle model from an input stream.
         * @param stream Input stream.
         * @param catalog Particle specifications catalog.
         * @return Coarse grained particle model.
         */
        static cg_ptr_t readFrom(std::istream& stream, 
                                 const spec_catalog_ptr_t& catalog);
        
        /**
         * Perform a task with protonatable beads. The given task must expose the 
         * operator
         * <code>
         *  R operator () (const std::vector<prot_bead_ptr_t>& protonatableBeads);
         * </code>
         * @param task Task. Maybe a lambda expression.
         * @return Result.
         */
        template <typename R, typename TASK>
        R doWithProtBeads(const TASK& task) { return task(discrete_, continouos_); }
        
    private:
        
        std::vector<dprot_bead_ptr_t> discrete_;
        std::vector<cprot_bead_ptr_t> continouos_;
    };    
    
    /**
     * Writes coarse grained particle model to output stream.
     * @param stream Output stream.
     * @param cg Coarse grained particle model.
     * @return Output stream.
     */
    std::ostream& operator << (std::ostream& stream, const CoarseGrained& cg);
    
}

#endif /* COARSE_GRAINED_HPP */

