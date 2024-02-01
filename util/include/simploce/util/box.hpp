/*
 * File:   box.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 29, 2019, 5:29 PM
 */

#ifndef BOX_HPP
#define BOX_HPP

#include "simploce/conf/u-conf.hpp"
#include "../types/u-types.hpp"
#include <cassert>
#include <array>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace simploce {
    
    /**
     * 3D Box with arbitrary box side lengths.
     * @param T value type.
     */
    template <typename V = real_t>
    class Box {
    public:
    
        /**
         * Default constructable. Box with zero edge lengths.
         */
        Box() : lengths_{0.0, 0.0, 0.0} {}
        
        /**
         * Box with equal side length.
         * @param v Side length.
         */
        explicit Box(V v) : lengths_{v, v, v} {}
    
        /**
         * Box with arbitrary edge lengths.
         * @param lx Dimension in x-direction.
         * @param ly Dimension in y-direction.
         * @param lz Dimension in z-direction.
         */
        Box(V lx, V ly, V lz) : lengths_{lx, ly, lz} {}

        virtual ~Box() = default;

        /**
         * Returns box dimension (edge length) in x-direction.
         * return Box size.
         */
        V lengthX() const { return lengths_[0]; }
    
        /**
         * Returns box dimension (edge length) in y-direction.
         * return Box size.
         */
        V lengthY() const { return lengths_[1]; }
    
        /**
         * Returns box dimension (edge length) in z-direction.
         * return Box size.
        */
        V lengthZ() const { return lengths_[2]; }
    
        /**
         * Returns box dimension (edge length) for given dimension.
         * @param k Dimension (k < 3, k = 0 -> x, k = 1 -> y, k = 2 -> z).
         * @return Length.
         */
        V operator [] (std::size_t k) const { return lengths_[k]; }
    
        /**
         * Returns box edge length for given dimension.
         * @param k Dimension (k < 3, k = 0 -> x, k = 1 -> y, k = 2 -> z).
         * @return Length.
         */
        V operator () (std::size_t k) const { return lengths_[k]; }
    
        /**
         * Returns size of the box. This is the largest edge or side length.
         * @return Length.
         */
        virtual V size() const;
    
        /**
         * Returns volume.
         * @return Volume.
         */
        virtual V volume() const { return lengths_[0] * lengths_[1] * lengths_[2]; }

        /**
         * Returns position at the center of the box.
         * @return Center.
         */
        position_t center() const;
    
    protected:
    
        /**
         * Sets edge lengths of the box.
         * @param lx Dimension in x-direction.
         * @param ly Dimension in y-direction.
         * @param lz Dimension in z-direction.
         */
        void lengths(V lx, V ly, V lz) { lengths_ = std::array<V,3>{lx, ly, lz}; }
    
    
    private:

        template <typename TT>
        friend std::istream& operator >> (std::istream& stream, Box<TT>& box);
    
        std::array<V, 3> lengths_;
    
    };
    
    template <typename V>
    V Box<V>::size() const 
    {
        V maxLength = std::max(lengths_[0], lengths_[1]);
        return std::max(maxLength, lengths_[2]);
    }

    template <typename V>
    position_t
    Box<V>::center() const {
        return position_t{0.5 * lengthX(), 0.5 * lengthY(), 0.5 * lengthZ()};
    }

    /**
     * Writes box to output stream.
     * @param stream Output stream.
     * @param box Box.
     * @return Output stream.
     */
    template <typename V>
    std::ostream& operator << (std::ostream &stream, const Box<V>& box)
    {
        const int precision = 7;
        const int width = 15;
        const char space = ' ';
    
        stream.precision(precision);
        stream << std::setw(width) << box[0];
        stream << space << std::setw(width) << box[1];
        stream << space << std::setw(width) << box[2];
        return stream;
    }
  
        /**
         * Reads box from input stream.
         * @param stream Input stream.
         * @param box Box.
         * @return Input stream.
         */
    template <typename V>
    std::istream& operator >> (std::istream& stream, Box<V>& box)
    {
        V lengthX, lengthY, lengthZ;
        stream >> lengthX >> lengthY >> lengthZ;
        box.lengths(lengthX, lengthY, lengthZ);
        return stream;
    }
    
    template <typename V>
    std::ostream& operator << (std::ostream& stream, Box<V>& box)
    {
        stream.setf(std::ios::scientific);
        stream.precision(conf::PRECISION);
        stream << std::setw(conf::REAL_WIDTH) << box[0]
               << std::setw(conf::REAL_WIDTH) << box[1]
               << std::setw(conf::REAL_WIDTH) << box[2];
        return stream;
    }
    
}

#endif /* BOX_HPP */

