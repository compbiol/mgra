#ifndef SYM_MULTIHASHMAP_H_
#define SYM_MULTIHASHMAP_H_

namespace utility { 

template<class item_class, class Hash = std::hash<item_class> >
struct sym_multi_hashmap: public std::unordered_multimap<item_class, item_class, Hash> {
  typedef std::unordered_multimap<item_class, item_class, Hash> multi_hashmap;
  typedef typename multi_hashmap::const_iterator const_iterator;

  using multi_hashmap::begin;
  using multi_hashmap::end;
  using multi_hashmap::find;
  using multi_hashmap::equal_range;
  
  sym_multi_hashmap() {
    m_card = 0;
  }

  bool defined(item_class const & x) const {
    return (multi_hashmap::count(x) != 0);
  }

  bool defined(item_class const & x, item_class const & y) const {
    std::pair<const_iterator, const_iterator> range = multi_hashmap::equal_range(x);
    for (auto it = range.first; it != range.second; ++it) { 
      if (it->second == y) { 
      	return true; 
      } 
    } 
    return false;
  }

  void insert(item_class const & x, item_class const & y) {
    multi_hashmap::insert(std::make_pair(x, y));
    if (x != y) { 
      multi_hashmap::insert(std::make_pair(y, x));
    } 
    ++m_card;
  }

  const_iterator find(item_class const & x, item_class const & y) { 
    std::pair<const_iterator, const_iterator> range = multi_hashmap::equal_range(x);
    for (auto it = range.first; it != range.second; ++it) { 
      if (it->second == y) { 
      	return it; 
      } 
    } 
    return multi_hashmap::end();
  } 

  void erase(item_class const & x, item_class const & y) {
    if (!defined(x)) {
      std::cerr << "sym_multi_hashmap::erase() error: unmapped pair (" << x << "," << y << ")" << std::endl;
      abort();
    }

    const_iterator rem = find(x, y);
    assert(rem != multi_hashmap::end()); 
    multi_hashmap::erase(rem);

    if (x != y) {
      rem = find(y, x);
      assert(rem != multi_hashmap::end()); 
      multi_hashmap::erase(rem);
    }
    --m_card;
  }

  const item_class& operator[] (item_class const & x) const { //FIXME deleted
    auto ix = multi_hashmap::find(x);
    if (ix == multi_hashmap::end()) {
      std::cerr << "sym_multi_hashmap::operator[] error: undefined element " << x << std::endl;
      abort();
    }
    return ix->second;
  }

  size_t size() const {
    return m_card;
  }

  bool empty() const {
    return (m_card == 0);
  }

  void clear() {
    multi_hashmap::clear();
    m_card = 0;
  }

private:
  size_t m_card;
};

} 

#endif
