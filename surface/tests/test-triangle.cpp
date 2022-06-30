/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/16/22.
 */

#include "simploce/surface/triangle.hpp"
#include "simploce/surface/vertex.hpp"
#include <iostream>

using namespace simploce;

int main() {
    auto v1 = Vertex::create(position_t{0.0, 0.0, 0.0}, normal_t{});
    auto v2 = Vertex::create(position_t{0.0, 0.0, 2.0}, normal_t{});
    auto v3 = Vertex::create(position_t{0.0, 2.0, 0.0}, normal_t{});

    auto triangle = Triangle::create(v1, v2, v3);
    std::cout << "Area: " << triangle->area() << std::endl;

    return 0;
}