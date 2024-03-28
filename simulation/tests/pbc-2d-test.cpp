//
// Created by ajuffer on 2/1/24.
//

#include "simploce/simulation/bc.hpp"
#include "simploce/simulation/s-factory.hpp"
#include "simploce/util/direction.hpp"

using namespace simploce;

void test() {
    std::cout << "pdb-test" << std::endl;

    auto box = factory::box(2.0);
    std::cout << "Box size: " << box->size() << std::endl;
    bc_ptr_t bc = factory::pbc_2d(box, Direction::X, Direction::Y);

    position_t r1{0,0,0};
    real_t dz = 0.1;
    std::size_t N = 20;

    std::cout.setf(std::ios::scientific);
    for (std::size_t k = 1; k != N; ++k) {
        auto x = real_t(k) * dz;
        position_t r2{x, 0, 0};
        length_t R = norm<real_t>(r1-r2);
        auto r12 = bc->apply(r1, r2);
        length_t Rbc = norm<real_t>(r12);
        std::cout << "k, R, Rbc: " << k << ' ' << R << ' ' << Rbc << std::endl;
    }
}


int main(int argc, char** argv) {
    test();
}
