/*
 * Author: André H. Juffer, Biocenter Oulu, Finland.
 * Created on on 5/17/22.
 */

#include "./types/srf-types.hpp"

namespace simploce {

    /**
     * Point where two or more curves, lines, or edges meet. A vertex has coordinates and associated unit
     * normal vector.
     * @see <a href="https://en.wikipedia.org/wiki/Vertex_(geometry)">Wikipedia</a>
     */
    class Vertex {
    public:

        /**
         * Returns vertex.
         * @param r Position.
         * @param un Unit normal.
         * @return Vertex.
         */
        static vertex_ptr_t create(const position_t& r, const normal_t& un);

        /**
         * Constructor
         * @param r Location.
         * @param un Unit normal.
         */
        Vertex(position_t r, normal_t un);

        /**
         * Returns vertex index.
         * @return Index.
         */
        int index() const;

        /**
         * Returns x-, y-,or z-coordinate of position.
         * @param k One of {0,1,2}.
         * @return Value.
         */
        real_t operator [] (int k) const;

        /**
         * Returns position.
         * @return Position.
         */
        position_t position() const;

        /**
         * Returns unit normal.
         * @return Unit normal.
         */
        normal_t normal() const;

    private:

        friend class polyhedron_generator;
        friend class Polyhedron;

        static int INDEX_;

        // Sets position.
        void position_(const position_t& r);

        // Sets unit normal vector.
        void normal_(const normal_t& normal);

        int index_;     // Internal identifier.
        position_t r_;  // Location.
        normal_t un_;   // Unit normal vector.
    };
}