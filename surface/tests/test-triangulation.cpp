/*
 * Author: André H. Juffer.
 * Created on 15/05/2022, 16:25.
 *
 * Copyright (c) 2022 Biocenter Oulu, University of Oulu, Finland. All rights reserved.
 */

#include "simploce/surface/dotted-surface-generator.hpp"
#include "simploce/surface/triangulation.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/surface/edge.hpp"
#include "simploce/surface/face.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/logger.hpp"
#include <iostream>
#include <vector>

using namespace simploce;

int main() {
    util::Logger::changeLogLevel(util::Logger::LOGTRACE);

    std::vector<position_t> r;
    std::vector<radius_t> radii;

    radii.emplace_back(radius_t{0.12});
    radii.emplace_back(radius_t{0.18});

    r.emplace_back(position_t{0.0, 0.0, 0.0});
    r.emplace_back(position_t{0.13, 0.0, 0.0});

    radius_t radius = 3.0;
    //auto result = dotted_surface_generator::spherical(radius, 500);
    auto result = dotted_surface_generator::general(r, radii);
    auto dots = result.first;
    std::clog << "Number of dots: " << dots.size() << std::endl;
    std::clog << "Area: " << result.second << std::endl;

    triangulation gen;
    auto surface = gen.spherical(radius, 960);

    surface->doWithAll<void>([] (const std::vector<vertex_ptr_t>& vertices,
                                 const std::vector<face_ptr_t>& faces,
                                 const std::vector<edge_ptr_t>& edges) {
        std::clog << "Number of vertices: " << vertices.size() << std::endl;
        std::clog << "Number of faces (triangles): " << faces.size() << std::endl;
        std::clog << "Number of edges: " << edges.size() << std::endl;
        for (auto& edge: edges) {
            std::clog << "Edge id: " << edge->id() << std::endl;
        }
    });
    std::clog << "Euler characteristic: " << surface->eulerCharacteristic() << std::endl;
    std::clog << "Map onto dotted surface..." << std::endl;
    polyhedron_generator::mapOnto(dots, surface);
    std::clog << "Done." << std::endl;
    std::clog << "Triangulated surface area: " << surface->area()() << std::endl;

    std::ofstream ostream;
    util::open_output_file(ostream, "/wrk3/tests/dotted-surface.dat");
    ostream << dots.size() << std::endl;
    for (auto& dot: dots ) {
        ostream << dot << std::endl;
    }
    ostream.close();
    util::open_output_file(ostream, "/wrk3/tests/edges-triangulated-surface.dat");
    surface->writeEdges(ostream);
    ostream.close();

    util::open_output_file(ostream, "/wrk3/tests/triangulated-surface.dat");
    ostream << *surface << std::endl;
    ostream.close();

    surface->doWithAll<void>([] (const std::vector<vertex_ptr_t>& vertices,
                                      const std::vector<face_ptr_t>& faces,
                                      const std::vector<edge_ptr_t>& edges) {
        for (const auto& face : faces) {
            auto indices = face->indices();
            for (auto index: indices) {
                std::cout << " " << index;
            }
            auto pair = face->center();
            std::cout << " " << norm<real_t>(pair.first) << " " << norm<real_t>(pair.second) << std::endl;
        }
    });
}
