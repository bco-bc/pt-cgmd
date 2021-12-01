/*
 * Author: André H. Juffer.
 * Created on 26/11/2021, 17:01.
 *
 * Copyright (c) 2021 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef SIMULATION_DISTANCE_HPP
#define SIMULATION_DISTANCE_HPP

#include "analyzer.hpp"
#include "a-types.hpp"
#include "../simulation/s-types.hpp"

namespace simploce {

    /**
     * Monitor distance between two specific particles.
     */
    class Distance : public analyzer {
    public:
        /**
         * Creates distance analyzer.
         * @param analysisParameters Analysis parameters.
         * @param bc Boundary condition.
         * @return Distance analyzer.
         */
        static d_ptr_t
        create(const a_param_ptr_t& analysisParameters, const bc_ptr_t& bc);

        /**
         * Constructor.
         * @param id_i Particle identifier.
         * @param id_j Particle identifier.
         * @param bc Boundary condition.
         */
        Distance(const a_param_ptr_t& analysisParameters, bc_ptr_t bc);

        void
        perform(const p_system_ptr_t& particleSystem) override;

        /**
         * Returns results.
         * @return Distance. First of each pair is the step index, the second is the distance.
         */
        std::vector<std::pair<int, dist_t>>
        results();

    private:

        id_t id_i_;
        id_t id_j_;
        bc_ptr_t bc_;

        std::vector<std::pair<int, dist_t>> distances_{};
    };
}

#endif //SIMULATION_DISTANCE_HPP
