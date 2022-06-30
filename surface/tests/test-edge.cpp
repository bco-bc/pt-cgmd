/*
 * Author: André H. Juffer.
 * Created on 15/05/2022, 17:19.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/surface/edge.hpp"
#include "simploce/surface/vertex.hpp"

using namespace simploce;

int main() {

    vertex_ptr_t v1 = Vertex::create(position_t{1.0, 0.0, 0.0}, normal_t{1.0, 0.0, 0.0});
    vertex_ptr_t v2 = Vertex::create(position_t{0.0, 0.0, 1.0}, normal_t{0.0, 0.0, 1.0});

    auto edge = Edge::create(v1, v2);

    vertex_ptr_t v3 = Vertex::create(position_t{0.0, 1.0, 0.0}, normal_t{0.0, 1.0, 0.0});

    Edge::create(v1, v2);
    Edge::create(v2, v3);

    // Should fail.
    std::clog << "Following should fail: Edge::create(v1, v1):" << std::endl;
    Edge::create(v1, v1);

    return 0;
}