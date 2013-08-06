#ifndef SYM_MULTIHASHMAP_H_
#define SYM_MULTIHASHMAP_H_

#include <functional>

template<class item_class, class Hash = std::hash<item_class> >
struct sym_multi_hashmap: public std::unordered_multimap<item_class, item_class, Hash> {
  typedef std::unordered_multimap<item_class, item_class, Hash> multi_hashmap;
  typedef typename multi_hashmap::const_iterator const_iterator;

  using multi_hashmap::end;
  using multi_hashmap::find; //REMOVE
  using multi_hashmap::count;
  using multi_hashmap::equal_range;
  
  sym_multi_hashmap() {
    card = 0;
  }

  bool defined(const item_class& x) const {
    return count(x) != 0;
  }

  void insert(const item_class& x, const item_class& y) {
    multi_hashmap::insert(std::make_pair(x, y));
    if (x != y) { 
      multi_hashmap::insert(std::make_pair(y, x));
    } 
    ++card;
  }

  const_iterator find(const item_class& x, const item_class& y) { 
    std::pair<const_iterator, const_iterator> range = multi_hashmap::equal_range(x);
    for (auto it = range.first; it != range.second; ++it) { 
      if (it->second == y) { 
	return it; 
      } 
    } 
    return multi_hashmap::end();
  } 

  void erase(const item_class& x, const item_class& y) {
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
    --card;
  }

  const item_class& operator[] (const item_class& x) const { //FIXME deleted
    auto ix = find(x);
    if (ix == end()) {
      std::cerr << "sym_multi_hashmap::operator[] error: undefined element " << x << std::endl;
      abort();
    }
    return ix->second;
  }

  size_t size() const {
    return card;
  }

  bool empty() const {
    return (card == 0);
  }

  void clear() {
    multi_hashmap::clear();
    card = 0;
  }
private:
  size_t card;
};

#endif
