/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   value_t.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 9, 2019, 1:25 PM
 */

#ifndef VALUE_T_HPP
#define VALUE_T_HPP

#include "uconf.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace simploce {
    
    /**
     * Holds a number of given type V. The template parameter D (for discriminator) 
     * allows different value types, such as
     * <spec>
     *     using energy_t = simploce::value_t<double, 1>;
     *     using temperature_t = simploce::value_t<double, 2>;
     *     using volume_t = simploce::value_t<double, 3>;
     * </spec>
     * Thus, energy_t and temperature_t represent truly different types. As a 
     * consequence,
     * <pre>
     *    energy_t E{123};
     *    temperature_t T{300};
     *    auto result = E + T;
     * </pre>
     * makes no sense whatsoever. But 
     * <pre>
     *    volume_t V{123};
     *    auto p = n * k * T / V;
     * </pre>
     * does.
     * @param V value type (double, float).
     * @param D Discriminator.
     */
    template <typename V, int D>
    class value_t {
    public:
        
        /**
         * Default constructible.
         */
        value_t() : v_{0} {}
        
        /**
         * Constructor
         * @param v Value.
         */
        value_t(const V& v) : v_{v} {}                
        
        /**
         * Returns value.
         * @return Value.
         */
        V operator () () const { return v_; }     
        
        /**
         * Addition.
         * @param v Value.
         * @return This value.
         */
        value_t operator += (const value_t& v);
        
        /**
         * Subtraction.
         * @param v Value.
         * @return This value.
         */
        value_t operator -= (const value_t& v);
        
    private:
        
        V v_;
            
    };
    
    template <typename V, int D>
    value_t<V,D> value_t<V,D>::operator += (const value_t<V,D>& v)
    {
        v_ += v();
        return *this;
    }
    
    template <typename V, int D>
    value_t<V,D> value_t<V,D>::operator -= (const value_t<V,D>& v)
    {
        v_ -= v();
        return *this;
    }

    /**
     * Addition
     * @param v1 Value
     * @param v1 Value
     * @return Result, v1 + v2
     */
    template <typename V, int D>
    value_t<V,D> operator + (const value_t<V,D>& v1, const value_t<V,D>& v2)
    {
        return value_t<V,D>{v1() + v2()};
    }
    
    /**
     * Subtraction.
     * @param v1 Value.
     * @param v1 Value.
     * @return Result, v1 - v2.
     */
    template <typename V, int D>
    value_t<V,D> operator - (const value_t<V,D>& v1, const value_t<V,D>& v2)
    {
        return value_t<V,D>{v1() - v2()};
    }
    
    /**
     * Multiplication by number from the right.
     * @param v Value
     * @param number Number.
     * @return Result.
     */
    template <typename V, int D, typename N>
    value_t<V,D> operator * (const value_t<V,D>& v, N number)
    {
        return value_t<V,D>{v() * number};
    }
    
    /**
     * Multiplication by number from the left.
     * @param v Value
     * @param number Number.
     * @return Result.
     */
    template <typename N, typename V, int D>
    value_t<V,D> operator * (const N& number, const value_t<V,D>& v)
    {
        return value_t<V,D>{v() * number};
    }
    
    /**
     * Multiplication of values of different types. Usage includes the following
     * value_
     * @param v1 Value
     * @param v2 Value
     * @return Result of type V.
     */
    template <typename V, int D1, int D2>
    V operator * (const value_t<V,D1>& v1, const value_t<V,D2>& v2)
    {
        return v1() * v2();
    }
    
    /**
     * Division of number by value.
     * @param number Number.
     * @param v Value
     * @return Result.
     */
    template <typename V, int D>
    value_t<V,D> operator / (const V& number, const value_t<V,D>& v)
    {
        return value_t<V,D>{number / v()};
    }
    
    /**
     * Division of value by number.
     * @param v Value
     * @param number Number
     * @return Result.
     */
    template <typename V, int D>
    value_t<V,D> operator / (const value_t<V,D>& v, const V& number)
    {
        return value_t<V,D>{v() / number};
    }
    
    /*
     * Division
     * @param v1 Value
     * @param v2 Value
     * @return Result of type V.
     */
    template <typename V, int D1, int D2>
    V operator / (const value_t<V,D1>& v1, const value_t<V,D2>& v2) {
        return v1() / v2();
    }
    
    /**
     * Equality operator. Returns true only if v1() and v2() are numerically equal.
     * @param v1 Value.
     * @param v2 Value.
     * @return Result.
     */
    template <typename V, int D>
    bool operator == (const value_t<V,D>& v1, const value_t<V,D>& v2)
    {
        V eps = std::numeric_limits<V>::epsilon();
        return (std::fabs(v1()-v2()) < eps);
    }
    
    /**
     * Comparison, smaller then or equal to.
     * @param v Value.
     * @param number Number.
     * @return Result of v <= number.
     */
    template <typename V, int D>
    bool operator <= (const value_t<V,D>& v, V number)
    {
        return v() <= number;
    }
    
    /**
     * Comparison, greater then or equal to
     * @param v Value.
     * @param number Number.
     * @return Result of v >= number.
     */
    template <typename V, int D>
    bool operator >= (const value_t<V,D>& v, V number)
    {
        return v() >= number;
    }
    
    /**
     * Writes value to an output stream, , using scientific 
     * notation with precision 7 and width 15
     * @param stream Output stream.
     * @param v Value.
     * @return Output stream.
     */
    template <typename V, int D>
    std::ostream& operator << (std::ostream& stream, const value_t<V,D>& v)
    {
        const int precision = conf::PRECISION;
        const int width = conf::WIDTH;

        stream.precision(precision);
        stream.setf(std::ios::scientific);
        stream << std::setw(width) << v();
        return stream;
    }
    
}

#endif /* VALUE_T_HPP */

