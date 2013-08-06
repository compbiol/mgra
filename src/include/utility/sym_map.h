/*
** Module: Symmetric maps ver. 1.3
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
  typedef std::map<item_class, item_class, order> map_t;
  typedef typename map_t::const_iterator const_iterator;

  using map_t::end;
  using map_t::find;

  sym_map() {
    card = 0;
  }

  bool defined(const item_class& x) const {
    return find(x) != end();
  }

  void insert(const item_class& x,const item_class& y) {
    map_t::insert(std::make_pair(x, y));
    map_t::insert(std::make_pair(y, x));
    ++card;
  }

  void erase(const item_class& x) {
    if (!defined(x)) {
      std::cerr << std::endl << "symmap::erase() error: unmapped element (" << x << ")" << std::endl;
      abort();
    }
    map_t::erase(map_t::find(x));
    map_t::erase(x);
    --card;
  }

  void erase(const item_class& x, const item_class& y) {
    if(!defined(x) || map_t::operator[](x) != y) {
      std::cerr << "symmap::erase() error: unmapped pair (" << x << "," << y << ")" << std::endl;
      abort();
    }
    map_t::erase(x);
    map_t::erase(y);
    --card;
  }

  size_t size() const {
    return card;
  }

  bool empty() const {
    return (card == 0);
  }

  void clear() {
    map_t::clear();
    card = 0;
  }

 private:
  size_t card;
};

#endif
