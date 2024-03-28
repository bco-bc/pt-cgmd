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
         * Returns subsets of a list of items. The number of subsets is decided by
         * the number of available hardware threads.
         * @param items Items
         * @return Subsets of items.
         */
        template <typename T, template <typename, typename... Args> class CONT>
        std::vector<std::vector<T>> 
        makeSubLists(const CONT<T>& items) {
            using sublists_t = std::vector<std::vector<T>>;

            static util::Logger logger{"util::makeSubLists"};
            logger.trace("Entering.");

            if ( items.empty() ) {
               util::logAndThrow(logger, "Empty list of items. Cannot create subsets.");
            }            
            
            // Subsets.
            sublists_t subSets{};
            
            std::size_t counter = 0;      
            static const std::size_t numberOfSubSets = std::thread::hardware_concurrency();
            logger.debug("Number of available threads: " + std::to_string(numberOfSubSets));
            std::size_t numberOfItemsPerSubList = items.size() / numberOfSubSets;
            for (std::size_t k = 0; k != numberOfSubSets; ++k) {
                std::vector<T> single{};  // A single subset of items.
                std::size_t n = 0;
                while (counter != items.size() && n != numberOfItemsPerSubList) {                    
                    single.push_back(items[counter]);
                    counter += 1;
                    n += 1;
                }
                subSets.push_back(single); // Store the subset.
            }
    
            // Add remaining items, if any, to the last subset.
            if ( counter < items.size() ) {
                auto& last = *(subSets.end() - 1);
                while ( counter != items.size() ) {
                    last.push_back(items[counter]);
                    counter += 1;
                }
            }
            logger.debug(std::to_string(subSets.size()) + ": Number of subsets.");
            
            // Done.
            logger.trace("Leaving.");
            return subSets;
        }

        /**
         * Return unit vector with random components.
         * @return Unit vector.
         */
        dist_vect_t randomUnit();

        /**
         * Integrates real-valued function from a to b using a simple numerical procedure.
         * @param FUNCTION A real-valued function with a real-valued argument. May be a lambda expression.
         * @param a Lower bound.
         * @param b Upper bound.
         * @return Result.
         *
        real_t
        integrate(real_t FUNCTION(real_t),
                  real_t a,
                  real_t b);
                  */

        /**
         * Integrates real-valued function from a to b using a simple numerical procedure.
         * @param FUNCTION A real-valued function with a real-valued argument. May be a lambda expression.
         * @param a Lower bound.
         * @param b Upper bound.
         * @return Result.
         */
        real_t
        integrate(const std::function<real_t (real_t x)>& function,
            real_t a,
            real_t b);

    }
}

#endif /* UTIL_HPP */

