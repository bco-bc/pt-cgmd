/*
 * File:   util.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 19, 2019, 5:23 PM
 */

#ifndef UTIL_HPP
#define UTIL_HPP

#include "simploce/types/u-types.hpp"
#include "simploce/util/logger.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <ctime>
#include <random>
#include <cstdlib>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfenv>
#include <future>
#include <vector>
#include <utility>
#include <cassert>
#include <stdexcept>

namespace simploce {
    namespace util {

        /**
         * Signum.
         * T is any type that supports constructor T{0} and the < operator.
         * @param val Some value of type T.
         * @return -1 (negative number), 1 (positive number), or 0 (zero).
         * @see <a href="http://en.wikipedia.org/wiki/Sign_function">Sign function</a>
         */
        template <typename T>
        int sgn(T val) {
            return (T{0} < val) - (val < T{0});
        }

        /**
         * Returns new identifier.
         * @return Identifier.
         */
        id_t generateId();

        /**
         * Convert to identifier type.
         * @param v
         * @return Identifier of type id_t
         */
        id_t toId(const std::string& v);
    
        /**
         * Returns seed value. Uses current time.
         * @return Seed value.
         */
        template <typename V>
        V seedValue() {
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
        int nint(real_t val);
            
        /**
         * Random number generator in a given range.
         * @param max Max value .
         * @return Random value between 0 and max, uniformly sampled.
         */
        real_t randomUniform(real_t min = 0, real_t max = 1.0);
        
         /**
         * Random number generator in the range [0, 1.0].
         * @return Random value.
         */
        real_t random();

       /**
         * Waits for all future results before returning.
         * @tparam T Result type stored in future.
         * @param futures Futures.
         * @return Results of futures.
         */
        template<typename T>
        std::vector<T> waitForAll(std::vector<std::future<T> >& futures) {
            std::vector<T> results;
            for (auto& f : futures)
                results.push_back(f.get());
            return results;
        }       
        
        /**
         * Returns string representation of t.
         * @param t Simple type.
         * @return String.
         */
        template <typename T>
        std::string to_string(const T& t) {
            return boost::lexical_cast<std::string, T>(t);
        }

        // std::string toString(const simploce::id_t& id);
        
        /**
         * Returns sublists of a list of items. The number of sublists is decided by
         * the number of available hardware threads.
         * @param items Items
         * @return Sub lists of items.
         */
        template <typename T, template <typename, typename... Args> class CONT>
        std::vector<std::vector<T>> 
        makeSubLists(const CONT<T>& items) {
            using sublists_t = std::vector<std::vector<T>>;
            
            static util::Logger logger{"util::makeSubLists"};
           if ( items.empty() ) {
               util::logAndThrow(logger, "Empty list of items. Cannot create sublists.");
            }            
            
            // Sublists.
            sublists_t subLists{};
            
            std::size_t counter = 0;      
            static const std::size_t nsublists = std::thread::hardware_concurrency();
            logger.debug("Number of available threads: " + std::to_string(nsublists));
            std::size_t numberOfItemsPerSubList = items.size() / nsublists;                                                              
            for (std::size_t k = 0; k != nsublists; ++k) {
                std::vector<T> single{};  // A single sublist of items.
                std::size_t n = 0;
                while (counter != items.size() && n != numberOfItemsPerSubList) {                    
                    single.push_back(items[counter]);
                    counter += 1;
                    n += 1;
                }
                subLists.push_back(single); // Store the sublist.               
            }
    
            // Add remaining items, if any, to the last sublist.
            if ( counter < items.size() ) {
                auto& last = *(subLists.end() - 1);
                while ( counter != items.size() ) {
                    last.push_back(items[counter]);
                    counter += 1;
                }
            }
            
            // Done.
            return subLists;
        }

        /**
         * Return unit vector with random components.
         * @return Unit vector.
         */
        dist_vect_t randomUnit();
    }
}

#endif /* UTIL_HPP */

