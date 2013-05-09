#ifndef SYM_MULTIHASHMAP_H_
#define SYM_MULTIHASHMAP_H_

#include <unordered_map>
#include <utility>
#include <functional>
#include <iostream>

template<class item_class, class Hash = std::hash<item_class> >
struct sym_multi_hashmap: public std::unordered_multimap<item_class, item_class, Hash> {
  typedef std::unordered_multimap<item_class, item_class, Hash> multi_hashmap;
  typedef typename multi_hashmap::const_iterator const_iterator;

  using multi_hashmap::end;
  using multi_hashmap::find;
  using multi_hashmap::equal_range;
  
  sym_multi_hashmap() {
    card = 0;
  }

  bool defined(const item_class& x) const {
    return find(x) != end();
  }

  void insert(const item_class& x, const item_class& y) {
    multi_hashmap::insert(std::make_pair(x, y));
    multi_hashmap::insert(std::make_pair(y, x));
    ++card;
  }

  void erase(const item_class& x) {
    multi_hashmap::erase(multi_hashmap::find(x));
    multi_hashmap::erase(x);
    --card;
  }
  

  const item_class& operator[] (const item_class& x) const {
    auto ix = find(x);
    if (ix == end()) {
      std::cerr << "sym_multi_hashmap::operator[] error: undefined element " << x << std::endl;
      abort();
    }
    return ix->second;
  }

	
  void erase(const item_class& x, const item_class& y) {
    if(!defined(x)) {
      std::cerr << "sym_multi_hashmap::erase() error: unmapped pair (" << x << "," << y << ")" << std::endl;
      abort();
    }
    multi_hashmap::erase(x);
    multi_hashmap::erase(y);
    --card;
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

  void pop(item_class& x, item_class& y) {
    if (!card) {
      std::cerr << "sym_multi_hashmap::extract error: empty sym_multi_hashmap!\n";
      abort();
    }
    x = multi_hashmap::begin()->first;
    y = multi_hashmap::begin()->second;
    multi_hashmap::erase(x);
    multi_hashmap::erase(y);
    --card;
  }

private:
  size_t card;
};

#endif
