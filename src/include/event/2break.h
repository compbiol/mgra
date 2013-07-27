/* 
** Module: 2-Breaks and Transformations support
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008 - 2013 by Max Alekseyev <maxal@cse.sc.edu>
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#ifndef TWOBREAK_H_
#define TWOBREAK_H_

#include <string>
#include <utility>

template<class mcolor_t>
struct TwoBreak {
  typedef std::pair<vertex_t, vertex_t> arc_t;
  typedef typename mcolor_t::citer citer; 
	
  TwoBreak(const arc_t& a1, const arc_t& a2, const mcolor_t& Q)
  : arcs({a1, a2}) 
  , mcolor(Q) 
  {
  }

  TwoBreak(const vertex_t& x1, const vertex_t& y1, const vertex_t& x2, const vertex_t& y2, const mcolor_t& Q)
  : arcs({arc_t(x1, y1), arc_t(x2, y2)})
  , mcolor(Q) 
  {
  }
    
  inline TwoBreak inverse() const {
    return TwoBreak(arcs[0].first, arcs[1].first, arcs[0].second, arcs[1].second, mcolor);
  }

  inline arc_t get_arc(size_t index) const { 
    assert(index < 2); 
    return arcs[index];
  } 

  inline mcolor_t get_mcolor() const { 
    return mcolor; 
  } 

  inline citer cbegin_mcolor() const { 
    return mcolor.cbegin();
  } 

  inline citer cend_mcolor() const { 
    return mcolor.cend();
  } 

  void normalize();
  void apply_single(partgraph_t& SG) const;

  friend std::ostream& operator << (std::ostream& os, const TwoBreak& t) {
    os << "(" << t.arcs[0].first << "," << t.arcs[0].second << ")x(" 
	<< t.arcs[1].first << "," << t.arcs[1].second << "):{" 
	<< genome_match::mcolor_to_name(t.mcolor) << "}";
    return os;
  }
 
private: 
  arc_t arcs[2]; // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
  mcolor_t mcolor; 
};

template<class mcolor_t>
void TwoBreak<mcolor_t>::apply_single(partgraph_t& SG) const {
  for(size_t i = 0; i < 2; ++i) {
    if (arcs[i].first != Infty || arcs[i].second != Infty) {
      SG.erase(arcs[i].first, arcs[i].second);
    }
  }
	
  if (arcs[0].first != Infty || arcs[1].first != Infty) {
    SG.insert(arcs[0].first, arcs[1].first);
  }

  if (arcs[0].second != Infty || arcs[1].second != Infty) {
    SG.insert(arcs[0].second, arcs[1].second);
  }
} 

template<class mcolor_t>
void TwoBreak<mcolor_t>::normalize() {
  while (true) {
    if (arcs[0].first > arcs[0].second) {
      arcs[0] = arc_t(arcs[0].second, arcs[0].first);
      arcs[1] = arc_t(arcs[1].second, arcs[1].first);
    } else if (arcs[0].first > arcs[1].first || arcs[0].first > arcs[1].second ) {
      arc_t temp = arcs[0];
      arcs[0] = arcs[1];
      arcs[1] = temp;
    } else { 
      break; 
    } 
  }
}

#endif

