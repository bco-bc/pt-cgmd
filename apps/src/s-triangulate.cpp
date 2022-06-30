/*
 * Triangulates the surface of a particle system.
 * Author: André H. Juffer, Biocenter Oulu, University of Oulu, Finland.
 *
 * Created on June 3, 2022.
 */

#include "simploce/simulation/s-factory.hpp"
#include "simploce/particle/particle-system.hpp"
#include "simploce/particle/particle.hpp"
#include "simploce/particle/particle-spec.hpp"
#include "simploce/surface/dotted-surface-generator.hpp"
#include "simploce/surface/triangulation.hpp"
#include "simploce/surface/polyhedron.hpp"
#include "simploce/util/logger.hpp"
#include "simploce/util/file.hpp"
#include "simploce/util/program-options.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace simploce;
namespace po = boost::program_options;

static void optimizeSurface(const polyhedron_ptr_t& polyhedron) {
    util::Logger logger("optimizeSurface");
    logger.trace("Entering.");
    logger.warn("Not yet implemented.");
    logger.trace("Leaving");
}

static std::pair<polyhedron_ptr_t, std::vector<position_t>>
triangulate(const p_system_ptr_t& particleSystem, std::size_t numberOfTriangles, bool optimize) {
    static util::Logger logger("triangulate()");
    logger.trace("Entering");

    using result_t = std::pair<polyhedron_ptr_t, std::vector<position_t>>;
    auto result = particleSystem->doWithAll<result_t>([numberOfTriangles] (std::vector<p_ptr_t>& particles) {
        std::vector<position_t> positions{};
        std::vector<radius_t> radii{};
        for (const auto& p: particles) {
            positions.emplace_back(p->position());
            radii.emplace_back(p->spec()->radius());
        }
        auto result = simploce::dotted_surface_generator::general(positions, radii);
        auto dots = result.first;
        triangulation triangulator;
        auto polyhedron = triangulator.spherical(radius_t{1.0}, numberOfTriangles);
        polyhedron_generator::mapOnto(dots, polyhedron);
        logger.info(std::to_string(result.second()) + ": Area dotted surface.");
        logger.info(std::to_string(polyhedron->area()()) + ": Area triangulated surface.");
        return std::make_pair(polyhedron, dots);
    });

    // allow optimization of surface?
    if (optimize) {
        optimizeSurface(result.first);
    }

    logger.trace("Leaving");
    return std::move(result);
}

int main(int argc, char *argv[]) {
    util::Logger logger("s-triangulate::main");
    util::Logger::changeLogLevel(util::Logger::LOGWARN);

    std::string fnTriangulatedSurface{"triangulated-surface.dat"};
    std::string fnEdges("edges-triangulated-surface.dat");
    std::string fnDottedSurface("dotted-surface.dat");
    std::size_t numberOfTriangles = 3840;
    bool optimize(false);

    // Default parameters.
    auto param = factory::simulationParameters();

    po::options_description usage("Usage");
    util::addStandardOptions(usage);
    usage.add_options() (
            "pdb",
            "Input file format content is PDB."
    )(
            "fn-output-triangulated-surface,o",
            po::value<std::string>(&fnTriangulatedSurface),
            "Output file name of triangulated surface. Default is 'triangulated-surface.dat'."
    )(
            "fn-edges-triangulated-surface,e",
            po::value<std::string>(&fnEdges),
            "Output file name of edges triangulated surface. Default is 'edges-triangulated-surface.dat'"
    )(
            "fn-dotted-surface",
            po::value<std::string>(&fnDottedSurface),
            "Output file name of dotted surface. Default is 'dotted-surface.dat'"
    )(
            "number-of-triangles",
            po::value<std::size_t>(&numberOfTriangles),
            "Request number of triangles (default is 3840)."
    )(
            "optimize",
            "Optimize surface, meaning that an attempt is done to make edges similar"
            " in length as much as possible. May be very slow."
    );

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, usage), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << "Triangulate surface of a macromolecule:" << std::endl;
        std::cout << usage << "\n";
        return 0;
    }
    util::verbose(vm);
    if (vm.count("optimize")) {
        optimize = true;
    }

    auto particleSystem = util::getParticleSystem(vm, param);
    logger.info(std::to_string(particleSystem->numberOfParticles()) + ": Number of particles.");

    // Triangulate.
    particleSystem->setOriginToCenterOfMass();
    auto pair = triangulate(particleSystem, numberOfTriangles, optimize);
    auto surface = pair.first;
    auto dots = pair.second;

    std::ofstream stream;
    util::open_output_file(stream, fnTriangulatedSurface);
    stream << *surface << std::endl;
    logger.info(fnTriangulatedSurface + ": Wrote triangulated surface.");
    stream.close();

    util::open_output_file(stream, fnEdges);
    surface->writeEdges(stream);
    logger.info(fnEdges + ": Wrote edges triangulated surface.");
    stream.close();

    util::open_output_file(stream, fnDottedSurface);
    stream << dots.size() << std::endl;
    for (auto& dot: dots ) {
        stream << dot << std::endl;
    }
    logger.info(fnDottedSurface + ": Wrote dotted surface.");
    stream.close();

    // Done.
    return EXIT_SUCCESS;
}
