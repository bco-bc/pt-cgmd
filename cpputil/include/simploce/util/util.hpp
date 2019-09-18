/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   util.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 19, 2019, 5:23 PM
 */

#ifndef UTIL_HPP
#define UTIL_HPP

#include <ctime>
#include <random>
#include <cstdlib>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfenv>
#include <future>
#include <vector>

namespace simploce {
    namespace util {
    
        /**
         * Returns seed value. Uses current time.
         * @return Seed value.
         */
        template <typename V>
        V seedValue()
        {
            std::time_t timer;
            V seed = time(&timer);
            srand(seed);
            seed += std::rand(); 
            return seed;
        }
        
        /**
         * Returns nearest integer.
         * V is the value type.
         * @param val Val.
         * @return Nearest integer.
         */
        template <typename V>
        int nint(V val)
        {
            std::fesetround(FE_TONEAREST);
            return std::nearbyint(val);
        }
            
        /**
         * Very -simple- random number generator.
         * V is the value type of real numbers.
         * @return Random value in the range [0,1].
         */
        template <typename V>
        V random()
        {
            static bool init = false;
            if ( !init ) {
                std::time_t timer;
                V seed = time(&timer);
                srand(seed);
                init = true;
            }
            int n = std::rand();
            V v = V(n) / RAND_MAX;
            return v;
        }
        
        /**
         * Waits for all future results before returning. T is result type stored 
         * in the future.
         * @param futures Futures.
         * @return Results of futures.
         */
        template<typename T>
        std::vector<T> waitForAll(std::vector<std::future<T> >& futures)
        {
            std::vector<T> results;
            for (auto& f : futures)
                results.push_back(f.get());
            return results;
        }            
            
    }
}

#endif /* UTIL_HPP */

