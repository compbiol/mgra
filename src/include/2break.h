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

#include "mpbgraph.h"

typedef std::pair<vertex_t, vertex_t> arc_t;

//template<class mcolor_t>
struct TwoBreak {
  TwoBreak() {
  };

  TwoBreak(const arc_t& a1, const arc_t& a2, const Mcolor& Q)
  : MultiColor(Q) {
    OldArc[0] = a1;
    OldArc[1] = a2;
  }

  TwoBreak(const vertex_t& x1, const vertex_t& y1, const vertex_t& x2, const vertex_t& y2, const Mcolor& Q)
  : MultiColor(Q) {
    OldArc[0] = std::make_pair(x1, y1);
    OldArc[1] = std::make_pair(x2, y2);
  }

  inline void revert(MBGraph& M) const {
    TwoBreak(OldArc[0].first, OldArc[1].first, OldArc[0].second, OldArc[1].second, MultiColor).apply(M);
  }

  inline void revertSingle(partgraph_t& SG) const {
    TwoBreak(OldArc[0].first, OldArc[1].first, OldArc[0].second, OldArc[1].second, MultiColor).applySingle(SG);
  }

  inline TwoBreak inverse() const {
    return TwoBreak(OldArc[0].first, OldArc[1].first, OldArc[0].second, OldArc[1].second, MultiColor);
  }

  bool islinear(MBGraph& M) const;

  bool apply(MBGraph& M, bool record = false) const;

  void applySingle(partgraph_t& SG) const;

  // make OldArc[0].first the smallest possible
  void normalize();

  friend ostream& operator << (ostream& os, const TwoBreak& t) {
    os << "(" << t.OldArc[0].first << "," << t.OldArc[0].second << ")x(" 
       << t.OldArc[1].first << "," << t.OldArc[1].second << "):{" << genome_match::mcolor_to_name(t.MultiColor) << "}";
    return os;
  }

  static std::list<TwoBreak> History;

  arc_t OldArc[2];     // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
  Mcolor MultiColor;
};

typedef std::list<TwoBreak> transform_t;

void ApplyAll(MBGraph& M, const transform_t& T, bool record = false);
void RevertAll(MBGraph& M, const transform_t& T);

#endif

