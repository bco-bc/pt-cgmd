#ifndef MAP2_HPP
#define MAP2_HPP

#include <map>
#include <vector>
#include <iostream>

namespace simploce {

  /**
   * A 2D map of values, similar to a matrix where the element value (of type V) is obtained from 
   * two keys values (of type K) instead of one key value as in regular maps. V must be default 
   * constructible.
   */
  template <typename K, typename V>
  class MatrixMap {
  public:

    /**
     * Default constructor.
     */
    MatrixMap();

    /**
     * Adds value. The order of keys matters.
     * @param key1 First key.
     * @param key2 Second key.
     */
    void add(K key1, K key2, V value);

    /**
     * Removes value. The order of keys matters. Removal will succeed only if the element actually
     * exists.
     * @param key1 First key.
     * @param key2 Second key.
     */
    void remove(K key1, K key2);

    /**
     * Returns a reference to the mapped value of the element with keys equivalent to key1 and 
     * key2. 
     * @param key1 First key.
     * @param key2 Second key.
     * @return Value.
     * @throws std::out_of_range if nonexistent.
     */
    V at(K key1, K key2) const;

    /**
     * Returns value.
     * @param keys Two key values.
     * @return Value.
     * @throws std::out_of_range if nonexistent.
     */
    V get(const std::pair<K,K>& keys) const;

    /**
     * Already in map?
     * @param keys Two key values. The order of keys matters.
     * @return Result.
     */
    bool contains(K key1, K key2) const;

    /**
     * Returns all key pairs in this matrix map.
     * @return Key pairs.
     */
    std::vector<std::pair<K, K>> keyPairs() const;

    /**
     * Completely clears this map's content.
     */
    void clear();

    /**
     * Is this map empty?
     * @return Result.
     */
    bool empty() const;

  private:

    using map_t = std::map<K,V>;    
    using map_map_t = std::map<K, map_t>;

    map_map_t cont_;  // "V=cont_(K1, K2)" corresponds to std::map<K1<std::map<K2,V>
  };

  template <typename K, typename V>
  MatrixMap<K,V>::MatrixMap() : cont_{}
  {
  }

  template <typename K, typename V>
  void MatrixMap<K,V>::add(K key1, K key2, V value)
  {
    auto iter = cont_.find(key1);
    if ( iter == cont_.end() ) {
      map_t rmap{};
      auto rpair = std::make_pair(key2, value);
      rmap.insert(rpair);
      auto cpair = std::make_pair(key1, rmap);
      cont_.insert(cpair);    
    } else {
      auto rpair = std::make_pair(key2, value);
      iter->second.insert(rpair);
    }
  }

  template <typename K, typename V>
  void MatrixMap<K,V>::remove(K key1, K key2)
  {
    auto iter = cont_.find(key1);
    if ( iter != cont_.end() ) {
      map_t& rmap = iter->second;
      rmap.erase(key2);
    }
  }
  
  template <typename K, typename V>
  V MatrixMap<K,V>::at(K key1, K key2) const
  {
    auto iter1 = cont_.find(key1);
    if ( iter1 == cont_.end() ) {
      throw std::out_of_range("No element associated with keys.");
    } else {
      const map_t& rmap = iter1->second;
      const auto& iter2 = rmap.find(key2);
      if ( iter2 == rmap.end() ) {
	throw std::out_of_range("No element associated with keys.");
      }
      return iter2->second;
    }
  }

  template <typename K, typename V>
  inline V MatrixMap<K,V>::get(const std::pair<K,K>& keys) const
  {
    return this->at(keys.first, keys.second); 
  }
  
  template <typename K, typename V>
  std::vector<std::pair<K, K>> MatrixMap<K,V>::keyPairs() const
  {
    std::vector<std::pair<K, K>> pairs;
    for (auto iter1 = cont_.begin(); iter1 != cont_.end(); ++iter1) {
      K key1 = iter1->first;
      const map_t& rmap = iter1->second;
      for (auto iter2 = rmap.begin(); iter2 != rmap.end(); ++iter2) {
	K key2 = iter2->first;
	auto pair = std::make_pair(key1, key2);
	pairs.push_back(pair);
      }
    }
    return pairs;
  }

  template <typename K, typename V>
  bool MatrixMap<K,V>::contains(K key1, K key2) const
  {
    auto iter1 = cont_.find(key1);
    if ( iter1 == cont_.end() ) {
      return false;
    } else {
      const map_t& rmap = iter1->second;
      const auto& iter2 = rmap.find(key2);
      return (iter2 == rmap.end() ? false : true);
    }
  }

  template <typename K, typename V>
  inline void MatrixMap<K,V>::clear()
  {
    cont_.clear();
  }

  template <typename K, typename V>
  inline bool MatrixMap<K,V>::empty() const
  {
    return cont_.empty();
  }

  /**
   * Writes map content to output stream.
   * @param stream Output stream.
   * @param ma Matrix map.
   * @return Output stream.
   */
  template <typename K, typename V>
  std::ostream& operator << (std::ostream& stream, const MatrixMap<K,V>& map)
  {
    if (!map.empty() ) {
      std::vector<std::pair<K, K>> pairs = map.keyPairs();
      for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {
	auto pair = *iter;
	V v = map.get(pair);
	stream << pair.first << " " << pair.second << " " << v;
	if (iter < pairs.end() - 1) {
	  stream << std::endl;
	}
      }
    }
    return stream;
  }

  /**
   * Writes map content to output stream.
   * @param stream Output stream.
   * @param ma Matrix map.
   * @return Output stream.
   */
  template <typename K>
  std::ostream& operator <<
  (std::ostream& stream, const MatrixMap<K,std::pair<double, double>>& map)
  {
    if (!map.empty() ) {
      std::vector<std::pair<K, K>> pairs = map.keyPairs();
      for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {
	auto pair = *iter;
	auto v = map.get(pair);
	stream << pair.first << " " << pair.second << " " << v.first << " " << v.second;
	if (iter < pairs.end() - 1) {
	  stream << std::endl;
	}
      }
    }
    return stream;
  }
  
}

#endif
