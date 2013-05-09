/*
** Module: Symmetric maps ver. 1.5
**
** This file is part of the
** Multiple Genome Rearrangements and Ancestors (MGRA)
** reconstruction software.
**
** Copyright (C) 2008,09 by Max Alekseyev <maxal@cse.sc.edu>
** 
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html
*/

#ifndef SYMMAP_H
#define SYMMAP_H

#include <map>
#include <utility>
#include <functional>
#include <iostream>

template<class item_class, class order = std::less<item_class> >
struct sym_map : public std::map<item_class, item_class, order> {
  typedef std::map<item_class, item_class, order> mymap;
  typedef typename mymap::const_iterator const_iterator;

  using mymap::end;
  using mymap::find;

  sym_map() {
    card = 0;
  }

  bool defined(const item_class& x) const {
    return find(x) != end();
  }

  void insert(const item_class& x,const item_class& y) {
    /*if (defined(x) || defined(y)) { //FIXME
      std::cerr << std::endl << "symmap::insert() error: redefining ";
      if (defined(x)) { 
	std::cerr << "(" << x << "," << mymap::operator[](x) << ") ";
      } 
      if (defined(y)) { 
	std::cerr << "(" << y << "," << mymap::operator[](y) << ") ";
      } 
      std::cerr << "with (" << x << "," << y << ")" << std::endl;
      abort();
    }*/
    mymap::insert(std::make_pair(x, y));
    mymap::insert(std::make_pair(y, x));
    ++card;
  }

  void erase(const item_class& x) {
    if (!defined(x)) {
      //std::cerr << std::endl << "symmap::erase() error: unmapped element (" << x << ")" << std::endl;
      abort();
    }
    mymap::erase(mymap::find(x));
    mymap::erase(x);
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
    if(!defined(x) || mymap::operator[](x) != y) {
      //std::cerr << "symmap::erase() error: unmapped pair (" << x << "," << y << ")" << std::endl;
      abort();
    }
    mymap::erase(x);
    mymap::erase(y);
    --card;
  }

  size_t size() const {
    return card;
  }

  bool empty() const {
    return (card == 0);
  }

  void clear() {
    mymap::clear();
    card = 0;
  }

  void pop(item_class& x, item_class& y) {
    if (!card) {
      std::cerr << "symmap::extract error: empty symmap!\n";
      abort();
    }
    x = mymap::begin()->first;
    y = mymap::begin()->second;
    mymap::erase(x);
    mymap::erase(y);
    --card;
  }

private:
  size_t card;
};

#endif
