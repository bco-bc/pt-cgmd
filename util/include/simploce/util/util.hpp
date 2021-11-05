/*
 * File:   util.hpp
 * Author: André H. Juffer, Biocenter Oulu.
 *
 * Created on August 19, 2019, 5:23 PM
 */

#ifndef UTIL_HPP
#define UTIL_HPP

#include "simploce/types/u-types.hpp"
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
         * Very -simple- random number generator. Not for production.
         * V is the value type of real numbers.
         * @return Random value in the range [0,1].
         */
        template <typename V>
        V random() {
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
        std::string toString(const T& t) {
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
            
            if ( items.empty() ) {
                throw std::domain_error("Empty list of items. Cannot create sublists.");
            }            
            
            // Sub lists.
            sublists_t subLists{};
            
            std::size_t counter = 0;      
            static const std::size_t nsublists = std::thread::hardware_concurrency();        
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

    }
}

#endif /* UTIL_HPP */

