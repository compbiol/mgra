/* 
** Module: 2-Breaks and Transformations support
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

#ifndef TWOBREAK_H
#define TWOBREAK_H

#include <string>
#include <utility>

#include "utility/sym_multi_hashmap.h"

typedef std::string vertex_t; //FIXME
typedef sym_multi_hashmap<vertex_t> partgraph_t; //FIXME
const vertex_t Infty = "oo"; //FIXME

template<class graph_t, class mcolor_t>
struct TwoBreak {
  typedef std::pair<vertex_t, vertex_t> arc_t;

  TwoBreak() {
  };

  TwoBreak(const arc_t& a1, const arc_t& a2, const mcolor_t& Q)
  : MultiColor(Q) {
    OldArc[0] = a1;
    OldArc[1] = a2;
  }

  TwoBreak(const vertex_t& x1, const vertex_t& y1, const vertex_t& x2, const vertex_t& y2, const mcolor_t& Q)
  : MultiColor(Q) {
    OldArc[0] = arc_t(x1, y1); 
    OldArc[1] = arc_t(x2, y2); 
  }
    
  inline void revert(graph_t& M) const {
    TwoBreak(OldArc[0].first, OldArc[1].first, OldArc[0].second, OldArc[1].second, MultiColor).apply(M);
  }

  inline void revert_single(partgraph_t& SG) const { 
    TwoBreak(OldArc[0].first, OldArc[1].first, OldArc[0].second, OldArc[1].second, MultiColor).apply_single(SG);
  }

  inline TwoBreak inverse() const {
    return TwoBreak(OldArc[0].first, OldArc[1].first, OldArc[0].second, OldArc[1].second, MultiColor);
  }

  bool is_linear(graph_t& M) const; 

  bool apply(graph_t& M, bool record = false) const;

  void apply_single(partgraph_t& SG) const;

  void normalize();

  friend std::ostream& operator << (std::ostream& os, const TwoBreak& t) {
    os << "(" << t.OldArc[0].first << "," << t.OldArc[0].second << ")x(" << t.OldArc[1].first << "," << t.OldArc[1].second << "):{" << genome_match::mcolor_to_name(t.MultiColor) << "}";
    return os;
  }
 
  //private: 
  arc_t OldArc[2];     // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
  mcolor_t MultiColor; //FIXME
};

// check if a 2-break creates a circular chromosome
template<class graph_t, class mcolor_t>
bool TwoBreak<graph_t, mcolor_t>::is_linear(graph_t& M) const {
  apply(M);

  for(int i = 0; i < 2; ++i) {
    vertex_t x;
    if (i == 0) {
      x = OldArc[0].first;
    } else { 
      x = OldArc[0].second;
    }

    if (x == Infty) {
      continue;
    }

    for(auto ic = MultiColor.cbegin(); ic != MultiColor.cend(); ++ic) {
      const partgraph_t& PG = M.get_local_graph(ic->first);
      bool circular = false;
      vertex_t y = M.get_adj_vertex(x);

      while (PG.defined(y)) { 
	y = PG[y];
	if (y != x) { 
	  y = M.get_adj_vertex(y);
	} 
	if (y == x) {
	  circular = true;
	  break;
	}
      }
		
      if (!circular && PG.defined(x)) {
	vertex_t y = x;
	while (PG.defined(y)) { 
	  y = PG[y];
	  if (y != x) {
	    y = M.get_adj_vertex(y);
	  } 
	  if (y == x) {
	    circular = true;
	    break;
	  }
	}
      }

      if (circular) {
	revert(M);
	return false;
      }
    }
  }

  revert(M);
  return true;
}

template<class graph_t, class mcolor_t>
bool TwoBreak<graph_t, mcolor_t>::apply(graph_t& M, bool record) const  {
  if (record) {
    M.insert_twobreak(*this);
  }
 
  for (auto ic = MultiColor.cbegin(); ic != MultiColor.cend(); ++ic) {
    for (size_t i = 0; i < 2; ++i) {
      if (OldArc[i].first != Infty && OldArc[i].second != Infty) {
	M.erase_edge(ic->first, OldArc[i].first, OldArc[i].second);
      } 
    }

    if(OldArc[0].first != Infty && OldArc[1].first != Infty) {
      M.add_edge(ic->first, OldArc[0].first, OldArc[1].first);
    }

    if(OldArc[0].second != Infty && OldArc[1].second != Infty) {
      M.add_edge(ic->first, OldArc[0].second, OldArc[1].second);
    }
  }
  return true;
}

template<class graph_t, class mcolor_t>
void TwoBreak<graph_t, mcolor_t>::apply_single(partgraph_t& SG) const {
  for(size_t i = 0; i < 2; ++i) {
    if (OldArc[i].first != Infty && OldArc[i].second != Infty) {
      SG.erase(OldArc[i].first, OldArc[i].second);
    }
  }

  if(OldArc[0].first != Infty && OldArc[1].first != Infty) {
    SG.insert(OldArc[0].first, OldArc[1].first);
  }

  if (OldArc[0].second != Infty && OldArc[1].second != Infty) {
    SG.insert(OldArc[0].second, OldArc[1].second);
  }
}

template<class graph_t, class mcolor_t>
void TwoBreak<graph_t, mcolor_t>::normalize() {
  while (true) {
    if (OldArc[0].first > OldArc[0].second) {
      OldArc[0] = arc_t(OldArc[0].second, OldArc[0].first);
      OldArc[1] = arc_t(OldArc[1].second, OldArc[1].first);
      continue;
    }

    if (OldArc[0].first > OldArc[1].first || OldArc[0].first > OldArc[1].second ) {
      arc_t temp = OldArc[0];
      OldArc[0] = OldArc[1];
      OldArc[1] = temp;
      continue;
    }

    break;
  }
}

#endif

