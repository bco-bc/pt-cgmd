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
        using lu_ptr_t = std::shared_ptr<lu_t>;

        /**
         * Constructor.
         * @param param BEM parameters.
         * @param nCol Number of collocation points.
         */
        BEMData(const param_ptr_t &param, std::size_t nCol);

        // Surface matrix.
        matrix_t S;

        // LU-decomposed S. Note that elements of S are replaced. Thus, no additional memory
        // is required.
        // see https://eigen.tuxfamily.org/dox/group__InplaceDecomposition.html
        lu_t lu;

        // Right-hand-side.
        vector_t b{};

        // Unknowns at collocation points.
        vector_t x{};

        // Dielectric constant of solute.
        real_t epsSolute{};

        // Dielectric constant of solvent.
        real_t epsSolvent{};

        // Ratio between the dielectric constants outside and inside the dielectric
        // boundary (surface).
        real_t epsRatio{};

        // Inverse Debye length.
        real_t ka{};
    };

}


#endif //BEM_BEM_DATA_HPP
