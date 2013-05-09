/* 
** Module: Equivalence relation ver. 1.3
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008,09 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#ifndef EQUIV_H
#define EQUIV_H

#include <map>
#include <functional>

template <class Item, class Cmp = std::less<Item> >
struct equivalence {
    typedef std::map<Item, Item, Cmp> map_t;

    // introduce a relation between two specified integers
    void addrel(const Item& x, const Item& y);

    void addrel(const std::pair<Item, Item>& p) { 
	addrel(p.first, p.second); 
    }

    // return `true' iff two integers are equivalent
    inline bool isequiv(const Item& x, const Item& y) { 
	return (operator[](x) == operator[](y));
    } 

    inline bool isequiv(const std::pair<Item, Item>& p) { 
	return isequiv(p.first, p.second); 
    }

    const Item& operator[](const Item& x);

    inline void insert(const Item& x) {
	operator[](x);
    }

    inline bool defined(const Item& x) const { 
	return (container.find(x) != container.end());
    } 

    void update();

    int classes();

    template<class eclass_t>
    void get_eclasses(std::map<Item, eclass_t, Cmp>& C);

private: 
    map_t container; 
};

template<class Item, class Cmp>
const Item& equivalence<Item,Cmp>::operator[] (const Item& x) {
    if (container.find(x) == container.end()) { 
	return container[x] = x; 
    } 

    Item y = x;
    while (container[y] != y) { 
	y = container[y];
    } 
    return (container[x] = y);
}

template<class Item, class Cmp>
void equivalence<Item,Cmp>::addrel(const Item& x, const Item& y) {
    Item z = operator[](x);
    Item t = operator[](y);

    if (container.key_comp()(z, t)) { 
	container[z] = t; 
    } else { 
	container[t] = z;
    } 
}

template<class Item, class Cmp>
void equivalence<Item,Cmp>::update() {
    for(auto mi = container.begin(); mi != container.end(); ++mi) { 
	operator[](mi->first);
    } 
}

template<class Item, class Cmp>
int equivalence<Item,Cmp>::classes() {
    int c = 0;
    for(auto mi = container.begin(); mi != container.end(); ++mi) { 
	if (operator[](mi->first) == mi->first) { 
	  ++c;
	} 
    } 
    return c;
}

template<class Item, class Cmp>
template<class eclass_t>
void equivalence<Item, Cmp>::get_eclasses(std::map<Item, eclass_t, Cmp>& C) {
    C.clear();
    for(auto mi = container.begin(); mi != container.end(); ++mi) { 
        C[operator[](mi->first)].insert(mi->first); 
   } 
}

/*
template<class Item, class Cmp>
template<class eclass_t>
void equivalence<Item,Cmp>::get_eclasses(map<Item,eclass_t,Cmp>& C) const {
    C.clear();
    for(typename map_t::const_iterator mi=begin();mi!=end();++mi)
        C[mi->second].insert(mi->first);
}
*/

#endif
