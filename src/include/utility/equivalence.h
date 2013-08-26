#ifndef EQUIV_H_
#define EQUIV_H_

/* 
** Equivalence relation.
**/

namespace utility { 

template <class Item, class Cmp = std::less<Item> >
struct equivalence {
  typedef std::map<Item, Item, Cmp> map_t;

  inline void addrel(const Item& x, const Item& y) {
    Item z = operator[](x);
    Item t = operator[](y);

    if (container.key_comp()(z, t)) { 
      container[z] = t; 
    } else { 
      container[t] = z;
    }
  } 

  inline void addrel(const std::pair<Item, Item>& p) { 
    addrel(p.first, p.second); 
  }

  // return `true' iff two integers are equivalent
  inline bool isequiv(const Item& x, const Item& y) { 
    return (operator[](x) == operator[](y));
  } 

  inline bool isequiv(const std::pair<Item, Item>& p) { 
    return isequiv(p.first, p.second); 
  }

  inline const Item& operator[](const Item& x) {
    if (container.count(x) == 0) { 
      return container[x] = x; 
    } 

    Item y = x;
    while (container[y] != y) { 
      y = container[y];
    } 
    return (container[x] = y);
  }

  inline void insert(const Item& x) {
    operator[](x);
  }

  inline bool defined(const Item& x) const { 
    return (container.find(x) != container.end());
  } 

  inline void update() {
    for (const auto & elem : container) { 
      operator[](elem.first);
    }
  } 

  size_t classes();

  template<class eclass_t>
  std::map<Item, eclass_t, Cmp> get_eclasses();

private: 
  map_t container; 
};

} 

template<class Item, class Cmp>
size_t utility::equivalence<Item, Cmp>::classes() {
  size_t count = 0;
  for(const auto& item : container) { 
    if (operator[](item.first) == item.first) { 
      ++count;
    } 
  } 
  return count;
}

template<class Item, class Cmp>
template<class eclass_t>
std::map<Item, eclass_t, Cmp> utility::equivalence<Item, Cmp>::get_eclasses() {
  std::map<Item, eclass_t, Cmp> classes;
  for(const auto& item : container) { 
    classes[operator[](item.first)].insert(item.first); 
  } 
  return classes;
}

#endif
