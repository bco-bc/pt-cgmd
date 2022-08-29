/*
 * Author: André H. Juffer.
 * Created on 26/05/2022, 22:57.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#ifndef BEM_BEM_DATA_HPP
#define BEM_BEM_DATA_HPP

#include "types/bem-types.hpp"
#include "Eigen/Dense"


namespace simploce {

    /**
     * BEM data holder.
     */
    struct BEMData {

        /**
         * Vector and matrix index type.
         */
        using index_t = Eigen::Index;

        /*
         * Type for square dense matrix.
         */
        using matrix_t = Eigen::MatrixXd;

        /**
         * Vector type.
         */
        using vector_t = Eigen::VectorXd;

        // Types for LU-decomposition.
        using lu_t = Eigen::PartialPivLU<Eigen::Ref<matrix_t>>;
        //using lu_ptr_t = std::shared_ptr<lu_t>;

        // Node
        using node_t = struct Node {
            Node() : r{}, nv{}, index{0}, area{0.0} {}
            position_t r;
            normal_t nv;
            index_t index;
            area_t area;
        };

        /**
         * Constructor.
         * @param param BEM parameters.
         * @param dimension Dimension of square surface matrix.
         */
        BEMData(const param_ptr_t &param, std::size_t dimension);

        // Nodes
        std::vector<node_t> nodes;

        // Surface matrix.
        matrix_t S;

        // LU-decomposed S. Note that elements of S are replaced. Thus, no additional memory
        // is required.
        // see https://eigen.tuxfamily.org/dox/group__InplaceDecomposition.html
        lu_t lu;

        // Right-hand-side.
        vector_t b{};

        // Unknowns.
        vector_t x{};

        // Dielectric constant of solute.
        real_t epsSolute{};

        // Dielectric constant of solvent.
        real_t epsSolvent{};

        // Ratio between the dielectric constants outside and inside the dielectric
        // boundary (surface).
        real_t epsRatio{};
    };

}


#endif //BEM_BEM_DATA_HPP
