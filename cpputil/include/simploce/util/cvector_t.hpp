/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cvector.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 8, 2019, 11:16 AM
 */

#ifndef CVECTOR_HPP
#define CVECTOR_HPP

#include "uconf.hpp"
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <array>

namespace simploce {
    
    /**
     * Cartesian vector with three components of given value type V. The
     * template parameter D (for discriminator) makes it possible to define various 
     * different Cartesian vector types, such as
     * <pre>
     *   using position_t = simploce::cvector<double, 1>;
     *   using velocity_t = simploce::cvector<double, 2>;
     * </pre>
     * Thus, position_t and velocity_t truly represent different types.
     * @see 
     * <a href="http://en.wikipedia.org/wiki/Cartesian_coordinate_system">
     *   Cartesian coordinate system
     * </a>
     * @see 
     * <a href="http://en.wikipedia.org/wiki/Euclidean_vector">
     *   Euclidian vector
     * </a>
     * @param V Value type (float or double, not complex).
     * @param D Discriminator.
     */
    template <typename V, int D>
    class cvector_t {
    public:    
    
        /**
         * Default constructible. Initialized with zero valued components.
         */
        cvector_t() : elems_{0.0, 0.0, 0.0} {}
        
        /**
         * Constructor.
         * @param v Vector. Its size must be at least 3.
         */
        cvector_t(const std::vector<V>& v) : elems_{v[0], v[1], v[2]} {}
            
        /**
         * Constructor.
         * @param x X-component.
         * @param y Y-component.
         * @param z Z-component.
         */
        cvector_t(V x, V y, V z) : elems_{x, y, z} {}
        
        cvector_t(const std::array<V,3>& a) : elems_{a} {}
    
        /**
         * Returns x-component.
         * @return Value.
         */
        V x() const { return elems_[0]; }
    
        /**
         * Returns y-component.
         * @return Value.
         */
        V y() const { return elems_[1]; }
    
        /**
         * Returns z-component.
         * @return Value.
         */
        V z() const { return elems_[2]; }
    
        /**
         * Returns modifiable value of vector component.
         * @param k Identifies component: k = 0,1,2 <-> x-, y-, z-component.
         * @return Value.
         */
        V& operator [] (std::size_t k) { return elems_[k]; }
    
        /**
         * Returns unmodifiable value of vector component.
         * @param k Identifies component: k = 0,1,2 <-> x-, y-, z-component.
         * @return Value.
         */
        V operator [] (std::size_t k) const { return elems_[k]; }
    
        /**
         * Addition.
         * @param v Cartesian vector.
         * @return This Cartesian vector.
         */
        cvector_t<V,D>& operator += (const cvector_t<V,D>& v);
    
        /**
         * Subtraction.
         * @param v Cartesian vector.
         * @return This Cartesian vector.
         */
        cvector_t<V,D>& operator -= (const cvector_t<V,D>& v);

        /**
         * Multiplication with a number.
         * @param t Number of type T.
         * @return This Cartesian vector.
         */
        template <typename T>
        cvector_t<V,D>& operator *= (T t); 
    
        /**
         * Division by a number,
         * @param t Number of type T. Must not be zero.
         * @return This Cartesian vector.
         */
        template <typename T>
        cvector_t<V,D>& operator /= (T t);     
    
        /**
         * Resets all components to zero.
         */
        void reset();
        
        /**
         * Returns all elements.
         * @return Elements.
         */
        std::array<V,3> toArray() { return elems_; }

    private:
        
        std::array<V, 3> elems_;
    };
    
    
    template <typename V, int D>
    cvector_t<V,D>& cvector_t<V,D>::operator += (const cvector_t<V,D>& v)
    {
        for (std::size_t k = 0; k != 3; ++k) {
            elems_[k] += v.elems_[k];
        }
        return *this;
    }
  
    template <typename V, int D>
    cvector_t<V,D>& cvector_t<V,D>::operator -= (const cvector_t<V,D>& v)
    {
        for (std::size_t k = 0; k != 3; ++k) {
            elems_[k] -= v.elems_[k];
        }
        return *this;
    }

    template <typename V, int D>
    template <typename T>
    cvector_t<V,D>& cvector_t<V,D>::operator *= (T t)
    {
        for (std::size_t k = 0; k != 3; ++k) {
            elems_[k] *= t;
        }
        return *this;
    }
  
    template <typename V, int D>
    template <typename T>
    cvector_t<V,D>& cvector_t<V,D>::operator /= (T t)
    {
        for (std::size_t k = 0; k != 3; ++k) {
            elems_[k] /= t;
        }
        return *this;
    }

    template <typename V, int D>
    void cvector_t<V,D>::reset()
    {
        for (std::size_t k = 0; k != 3; ++k) {
            elems_[k] = 0;
        }
    }
  
    /**
     * Adds two vectors of the -same- type.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @return Result: v1 + v2
     */
    template <typename V, int D>
    cvector_t<V,D> operator + (const cvector_t<V,D>& v1, const cvector_t<V,D>& v2)
    {
        cvector_t<V,D> v{};
        for ( std::size_t k = 0; k != 3; ++k ) {
            v[k] = v1[k] + v2[k];
        }
        return v;
    }
  
    /**
     * Subtracts two vectors of the -same- type.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @return Result: v1 - v2
     */
    template <typename V, int D>
    cvector_t<V,D> operator - (const cvector_t<V,D>& v1, const cvector_t<V,D>& v2)
    {
        cvector_t<V,D> v{};
        for ( std::size_t k = 0; k != 3; ++k ) {
            v[k] = v1[k] - v2[k];
        }
        return v;
    }

    /**
     * Multiplies vector with number from the left.
     * @param number Number of type T.
     * @param v Vector.
     * @return Result: number * v.
     */
    template <typename T, typename V, int D>
    cvector_t<V,D> operator * (T number, const cvector_t<V,D>& v)
    {
        cvector_t<V,D> result;
        for (std::size_t k = 0; k != 3; ++k) {
            result[k] = number * v[k];
        }
        return result;
    }
  
    /**
     * Multiplies vector with number from the right.
     * @param v Vector.
     * @param number Number of type T.
     * @return Result: v * number.
     */
    template <typename T, typename V, int D>
    inline cvector_t<V,D> operator * (const cvector_t<V,D>& v, T number)
    {
        return number * v;
    }
  
    /**
     * Division by number.
     * @param number Number of type T. Must not be zero.
     * @param v Vector.
     * @return Result: v / number.
     */
    template <typename T, typename V, int D>
    inline cvector_t<V,D> operator / (const cvector_t<V,D>& v, T number)
    {
        return ( T{1.0}/number ) * v;
    }
    
    /**
     * Equality operator. Returns true if (1) each component of v1 is numerically equal to 
     * the corresponding component of v2, and (2) v1 and v2 must be of the -same- type.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @return Result: v1 == v2.
     */
    template <typename V, int D>
    bool operator == (const cvector_t<V,D>& v1, const cvector_t<V,D>& v2)
    {
        V eps = std::numeric_limits<V>::epsilon();
        return ( std::fabs(v1[0] - v2[0]) < eps &&
                 std::fabs(v1[1] - v2[1]) < eps &&
                 std::fabs(v1[2] - v2[2]) < eps );
    }
  
    /**
     * Not equal operator. Returns true (1) if any component of v1 is numerically not 
     * equal to the corresponding component of v2, or (2) v1 and v2 represent different 
     * types.
     * @param v1 Vector 1.
     * @param v2 Vector 2.
     * @return Result: v1 != v2.
     */
    template <typename V, int D>
    inline bool operator != (const cvector_t<V,D>& v1, const cvector_t<V,D>& v2)
    {
        return ( v1 == v2 ? false : true );
    }
  
    /**
     * Returns norm (magnitude, length) of Cartesian vector.
     * @param v Cartesian vector.
     * @return Norm, a value of type T
     */
    template <typename T, typename V, int D>
    T norm(const cvector_t<V,D>& v)
    {
        return std::sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
    }

    /**
     * Returns square of norm (magnitude, length) of Cartesian vector.
     * @param v Cartesian vector.
     * @return Square norm, a value of type T.
     */
    template <typename T, typename V, int D>
    T norm2(const cvector_t<V,D>& v)
    {
        return v.x() * v.x() + v.y() * v.y() + v.z() * v.z(); 
    }
  
    /**
     * Inner (or dot) product between two Cartesian vectors, possibly of different types.
     * @param v1 Vector of type defined by V1.
     * @param v2 Vector of type defined by V2.
     * @return Result. A number of type T.
     * @see <a href="https://en.wikipedia.org/wiki/Dot_product">Dot product</a>
     */
    template <typename T, typename V1, int D1, typename V2, int D2>
    T inner(const cvector_t<V1,D1>& v1, const cvector_t<V2,D2>& v2)
    {
        T result{0};
        for (std::size_t k = 0; k != 3; ++k) {
            result += v1[k] * v2[k];
        }
        return result;      
    }
  
    /**
     * Returns angle between two Cartesian vectors, possibly of different types.
     * Use it as 
     * <pre>
     *     double a = simploce::angle<double>(v1, v2);
     * </pre>
     * @param v1 Vector of type defined by V1 and D1. Should have a non-zero norm.
     * @param v2 Vector of type defined by V2 and D2. Should have a non-zero norm.
     * @return Angle in rad (radian) in the range [0;pi].
     */
    template <typename T, typename V1, int D1, typename V2, int D2>
    T angle(const cvector_t<V1,D1>& v1, const cvector_t<V2,D2>& v2)
    { 
        T ip = inner<T>(v1, v2);
        T norm1 = norm<T>(v1);
        T norm2 = norm<T>(v2);
        T cos_angle = ip / (norm1 * norm2);
        return std::acos(cos_angle);
    }

    /**
     * Returns cross product between two Cartesion vectors, possibly of different types.
     * Usage:
     * <pre>
     *  simploce::cvector<double,2> cp = simploce::cross<double,2>(a, b);
     * </pre>
     * @param a Vector of type defined by V and D. Should have a non-zero norm.
     * @param b Vector of type defined by V and D. Should have a non-zero norm.
     * @return Cartesian vector v of type defined by V and D, v = a x b.
     */
    template <typename V, int D, int D1, int D2>
    cvector_t<V,D> cross(const cvector_t<V,D1>& a, const cvector_t<V,D2>& b)
    {
        V v0 = a[1] * b[2] - a[2] * b[1];
        V v1 = a[2] * b[0] - a[0] * b[2];
        V v2 = a[0] * b[1] - a[1] * b[0];
        return cvector_t<V,D>{v0,v1,v2};
    }

    /**
     * Writes each vector component to an output stream, using scientific 
     * notation with precision 7 and width 15.
     * @param stream Output stream.
     * @param v Cartesian vector.
     * @return Output stream.
     */
    template <typename V, int D>
    std::ostream& operator << (std::ostream& stream, const cvector_t<V,D>& v)
    {
        const int precision = conf::PRECISION;
        const int width = conf::WIDTH;
        const char space = conf::SPACE;
        
        stream.precision(precision);
        stream.setf(std::ios::scientific);
        stream << std::setw(width) << v.x();
        stream << space << std::setw(width) << v.y();
        stream << space << std::setw(width) << v.z();
        return stream;
    }
    
}

#endif /* CVECTOR_HPP */

