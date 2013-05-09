#ifndef SYM_HASHMAP_H_
#define SYM_HASHMAP_H_

#include <unordered_map>
#include <utility>
#include <functional>
#include <iostream>

template<class item_class, class Hash = std::hash<item_class> >
struct sym_hashmap: public std::unordered_map<item_class, item_class, Hash> {
  typedef std::unordered_map<item_class, item_class, Hash> hashmap;
  typedef typename hashmap::const_iterator const_iterator;

  using hashmap::end;
  using hashmap::find;

  sym_hashmap() {
    card = 0;
  }

  bool defined(const item_class& x) const {
    return find(x) != end();
  }

  void insert(const item_class& x, const item_class& y) {
    /*if (defined(x) || defined(y)) { //FIXME
      std::cerr << std::endl << "sym_hashmap::insert() error: redefining ";
      if (defined(x)) { 
	std::cerr << "(" << x << "," << myhash::operator[](x) << ") ";
      } 
      if (defined(y)) { 
	std::cerr << "(" << y << "," << myhash::operator[](y) << ") ";
      } 
      std::cerr << "with (" << x << "," << y << ")" << std::endl;
      abort();
    }*/
    hashmap::insert(std::make_pair(x, y));
    hashmap::insert(std::make_pair(y, x));
    ++card;
  }

  void erase(const item_class& x) {
    if (!defined(x)) {
      //std::cerr << std::endl << "symmap::erase() error: unmapped element (" << x << ")" << std::endl;
      abort();
    }
    hashmap::erase(hashmap::find(x));
    hashmap::erase(x);
    --card;
  }

  const item_class& operator [] (const item_class& x) const {
    auto ix = find(x);
    if (ix == end()) {
      //std::cerr << "symmap::operator[] error: undefined element " << x << std::endl;
      abort();
    }
    return ix->second;
  }

	
  void erase(const item_class& x, const item_class& y) {
    if(!defined(x) || hashmap::operator[](x) != y) {
      //std::cerr << "symmap::erase() error: unmapped pair (" << x << "," << y << ")" << std::endl;
      abort();
    }
    hashmap::erase(x);
    hashmap::erase(y);
    --card;
  }

  size_t size() const {
    return card;
  }

  bool empty() const {
    return (card == 0);
  }

  void clear() {
    hashmap::clear();
    card = 0;
  }

  void pop(item_class& x, item_class& y) {
    if (!card) {
      std::cerr << "sym_hashmap::extract error: empty symmap!\n";
      abort();
    }
    x = hashmap::begin()->first;
    y = hashmap::begin()->second;
    hashmap::erase(x);
    hashmap::erase(y);
    --card;
  }

private:
  size_t card;
};


#endif
