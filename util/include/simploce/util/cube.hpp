/*
 * File:   cube.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 29, 2019, 5:42 PM
 */

#ifndef CUBE_HPP
#define CUBE_HPP

#include "box.hpp"

namespace simploce {

    /**
     * A three-dimensional solid object bounded by six square faces, 
     * facets or sides, with three meeting at each vertex. The sides have equal length.
     * @param V value type.
     */
    template <typename V>
    class Cube : public Box<V> {
    public:
    
        /**
         * Constructs cube with zero edge length.
         */
        Cube() : Box<V>{} {}

        /**
         * Constructs box size with given edge or side length.
         * @param edgeLength Edge length.
         */
        explicit Cube(V v) : Box<V>{v} {}
    
        virtual ~Cube() = default;
        
        V size() const override { return this->lengthX(); }
    
        /**
         * Returns edge or side length.
         * @return Edge length.
         */
        V edgeLength() const { return this->lengthX(); }
    
    private:

        template <typename VV>
        friend std::istream & operator >> (std::istream& stream, Cube<VV>& cube);
    
        /**
         * Sets edge length.
         * @param length Edge length.
         */
        void edgeLength_(V edgeLength) { this->lengths(edgeLength, edgeLength, edgeLength); }
    
    };

    /**
     * Reads cube from input stream.
     * @param stream input stream. Should hold the cube's edge length.
     * @param cube Cube
     * @return Input stream.
     */
    template <typename V>
    std::istream& operator >> (std::istream& stream, Cube<V>& cube)
    {
        Box<V>& box = cube;
        stream >> box;
        return stream;
    }
  
    /**
     * Writes cube to output stream. Writes its edge length.
     * @param stream Output stream.
     * @param cube Cube
     * @return Output stream.
     */
    template <typename V>
    std::ostream& operator << (std::ostream& stream, const Cube<V>& cube)
    {
        const Box<V> box = cube;
        stream << box;
        return stream;
    }
    
}


#endif /* CUBE_HPP */

