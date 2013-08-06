/* 
** Module: Insertion/deletions operation support
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

#ifndef INSDEL_H_
#define INSDEL_H_

template<class mcolor_t>
struct InsDel {
  typedef std::pair<vertex_t, vertex_t> arc_t;
  typedef typename mcolor_t::citer citer; 
	
  InsDel(const arc_t& e, const mcolor_t& Q, bool is_del)
  : edge(e) 
  , mcolor(Q) 
  , is_deletion(is_del)
  {
  }

  InsDel(const vertex_t& x1, const vertex_t& y1, const mcolor_t& Q, bool is_del)
  : edge(arc_t(x1, y1)) 
  , mcolor(Q) 
  , is_deletion(is_del)
  {
  }
    
  inline InsDel inverse() const {
    return InsDel(edge, mcolor, !is_deletion);
  }

  inline bool is_deletion_oper() const { 
    return is_deletion;
  }
 
  inline arc_t get_edge() const { 
    return edge;
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

  friend std::ostream& operator << (std::ostream& os, const InsDel& t) {
    if (t.is_deletion) { 
	os << "Deletion "; 
    } else { 
	os << "Insertion ";
    } 
    os << "(" << t.edge.first << "," << t.edge.second << "):{" 
      << genome_match::mcolor_to_name(t.mcolor) << "}";    
    return os;
  }
 
private: 
  arc_t edge; // (x1, y1)
  mcolor_t mcolor; 
  bool is_deletion; 
};
#endif

