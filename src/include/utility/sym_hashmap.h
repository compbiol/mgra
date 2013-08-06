#ifndef SYM_HASHMAP_H_
#define SYM_HASHMAP_H_

#include <functional>

template<class item_class, class Hash = std::hash<item_class> >
  struct sym_hashmap: public std::unordered_map<item_class, item_class, Hash> {
  typedef std::unordered_map<item_class, item_class, Hash> hashmap;
  typedef typename hashmap::const_iterator const_iterator;

  using hashmap::end;
  using hashmap::find;

  sym_hashmap() {
    //card = 0;
  }

  bool defined(const item_class& x) const {
    return find(x) != end();
  }

  void insert(const item_class& x, const item_class& y) {
    hashmap::insert(std::make_pair(x, y));
    hashmap::insert(std::make_pair(y, x));
    //++card;
  }

  void erase(const item_class& x) {
    if (!defined(x)) {
      std::cerr << std::endl << "hashsymmap::erase() error: unmapped element (" << x << ")" << std::endl;
      abort();
    }
    hashmap::erase(hashmap::find(x));
    hashmap::erase(x);
    //--card;
  }
	
  /*size_t size() const {
    return card;
    }

    bool empty() const {
    return (card == 0);
    }
    private:
    size_t card;*/
};


#endif
