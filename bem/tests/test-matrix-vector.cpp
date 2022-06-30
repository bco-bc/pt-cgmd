/*
 * Author: André H. Juffer.
 * Created on 23/05/2022, 21:35.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/bem/types/bem-types.hpp"
#include "Eigen/Dense"

using matrix_t = Eigen::MatrixXd;

using namespace simploce;

int main() {
    std::size_t N = 3;
    matrix_t A(N, N);

    // Inverse.
    // Kreyszig p. 366, example 1.
    A(0, 0) = -1;
    A(0, 1) = 1;
    A(0, 2) = 2;
    A(1, 0) = 3;
    A(1, 1) = -1;
    A(1, 2) = 1;
    A(2, 0) = -1;
    A(2, 1) = 3;
    A(2, 2) = 4;
    std::cout << "A for inverse calculation: " << std::endl << A << std::endl;
    std::cout << "Inverse of A: " << std::endl << A.inverse() << std::endl;
}

