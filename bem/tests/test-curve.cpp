//
// Created by ajuffer on 7/15/22.
//

#include "simploce/bem/curve.hpp"
#include "simploce/surface/vertex.hpp"
#include <cstdlib>

using namespace simploce;

int main() {

    auto start = Vertex::create(position_t{0,0,1.0}, normal_t{0, 0, 1.0});
    auto end = Vertex::create(position_t{0, 1.0, 0.0}, normal_t{0, 1.0, 0.0});

    Curve curve(start, end);

    int N = 10;
    real_t dt = 0.1;
    for (auto k = 0; k != N + 1; ++k) {
        real_t t = k * dt;
        auto location = curve.point(t);
        std::cout << t << " " << location.first << location.second << std::endl;
    }

    return EXIT_SUCCESS;
}
