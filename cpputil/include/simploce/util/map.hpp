/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   map.hpp
 * Author: juffer
 *
 * Created on 22 September 2019, 14:26
 */

#ifndef MAP_HPP
#define MAP_HPP

#include <map>
#include <iostream>

namespace simploce {
    
    /**
     * Write map to an output stream.
     * @param stream Output stream.
     * @param map Map.
     * @return Output stream.
     */
    template <typename K, typename V>
    std::ostream& operator << (std::ostream& stream, std::map<K,V>& map)
    {
        for (auto iter = map.begin(); iter != map.end(); ++iter) {
            const auto& pair = *iter;
            stream << "(" << pair.first << ", " << pair.second << ")";
            if ( iter != --map.end()) {
                stream << std::endl;
            }
        }
        return stream;
    }
}

#endif /* MAP_HPP */

