#ifndef EQUIV_H_
#define EQUIV_H_

/* 
** Equivalence relation.
**/

namespace utility { 

template <class Item, class Cmp = std::less<Item> >
struct equivalence {
  typedef std::map<Item, Item, Cmp> map_t;

  inline void addrel(Item const & x, Item const & y) {
    Item z = operator[](x);
    Item t = operator[](y);

    if (container.key_comp()(z, t)) { 
      container[z] = t; 
    } else { 
      container[t] = z;
    }
  } 

  inline void addrel(std::pair<Item, Item> const & p) { 
    addrel(p.first, p.second); 
  }

  // return `true' iff two integers are equivalent
  inline bool isequiv(Item const & x, Item const & y) { 
    return (operator[](x) == operator[](y));
  } 

  inline bool isequiv(std::pair<Item, Item> const & p) { 
    return isequiv(p.first, p.second); 
  }

  inline Item const & operator[](Item const & x) {
    if (container.count(x) == 0) { 
      return container[x] = x; 
    } 

    Item y = x;
    while (container[y] != y) { 
      y = container[y];
    } 
    return (container[x] = y);
  }

  inline void insert(Item const & x) {
    operator[](x);
  }

  inline bool defined(Item const & x) const { 
    return (container.find(x) != container.end());
  } 

  inline void update() {
    for (auto const & elem : container) { 
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
  for(auto const & item : container) { 
    if (operator[](item.first) == item.first) { 
      ++count;
    } 
  } 
  return count;
}

template<class Item, class Cmp>
template<class eclass_t>
std::map<Item, eclass_t, Cmp> utility::equivalence<Item, Cmp>::get_eclasses() {
  std::map<Item, eclass_t , Cmp> classes;
  for(auto const & item : container) { 
    classes[operator[](item.first)].insert(item.first); 
  } 
  return classes;
}

#endif
