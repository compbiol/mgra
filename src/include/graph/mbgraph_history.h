/* 
** Module: Mutiple Breakpoint Graph with vec-T-consistent and T-consistent multicolors and save all events in evolution support
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

#ifndef MBGRAPH_HISTORY_H_
#define MBGRAPH_HISTORY_H_

#include "mbgraph_colors.h"
#include "2break.h"
#include "Insdel.h"
#include "TandemDuplication.h"

#define member(S,x) ((S).find(x)!=(S).end())

template<class mcolor_t>
struct mbgraph_with_history : public mbgraph_with_colors<mcolor_t> { 

	mbgraph_with_history(const std::vector<Genome>& genomes, const ProblemInstance& cfg)
	: mbgraph_with_colors<mcolor_t>(genomes, cfg) 
	{ 
	} 

	//2-break operations
	void apply_two_break(const TwoBreak<mcolor_t>& break2, bool record = false);
	void apply_single_two_break(size_t index, const TwoBreak<mcolor_t>& break2);  
	void apply_single_two_break(const TwoBreak<mcolor_t>& break2, partgraph_t& SG); //FIXME DEL 
	//bool is_linear(graph_t& M) const; 

	//Insertion/Deletion operations
	void apply_ins_del(const InsDel<mcolor_t>& insdel, bool record = true);
	void apply_single_ins_del(size_t index, const InsDel<mcolor_t>& insdel);  

	//(Reverse) tandem duplication operations
	void apply_tandem_duplication(const TandemDuplication<mcolor_t>& dupl, bool record = true);

	inline typename std::list<TwoBreak<mcolor_t> > get_history() const { //FIXME: DEL
		return break2_history; 
	} 

	inline typename std::list<TwoBreak<mcolor_t> >::const_iterator begin_history() const { 
		return break2_history.cbegin();
	} 
	
	inline typename std::list<TwoBreak<mcolor_t> >::const_iterator end_history() const { 
		return break2_history.cend(); 
	} 
private: 
	std::list<TwoBreak<mcolor_t> > break2_history;
	std::list<InsDel<mcolor_t> > insdel_history;		
	std::list<TandemDuplication<mcolor_t> > tandem_dupl_history;		
};



typedef std::list<TwoBreak<Mcolor> > transform_t;

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_two_break(const TwoBreak<mcolor_t>& break2, bool record) { 
  if (record) {
    break2_history.push_back(break2);
  }
 
  for (auto ic = break2.cbegin_mcolor(); ic != break2.cend_mcolor(); ++ic) {
    for (size_t i = 0; i < 2; ++i) {
      if (break2.get_arc(i).first != Infty || break2.get_arc(i).second != Infty) {
	erase_edge(ic->first, break2.get_arc(i).first, break2.get_arc(i).second);
      } 
    }

    if (break2.get_arc(0).first != Infty || break2.get_arc(1).first != Infty) {
      add_edge(ic->first, break2.get_arc(0).first, break2.get_arc(1).first);
    }

    if (break2.get_arc(0).second != Infty || break2.get_arc(1).second != Infty) {
      add_edge(ic->first, break2.get_arc(0).second, break2.get_arc(1).second);
    }
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_single_two_break(size_t index, const TwoBreak<mcolor_t>& break2) {
  for(size_t i = 0; i < 2; ++i) {
    if (break2.get_arc(i).first != Infty || break2.get_arc(i).second != Infty) {
      erase_edge(index, break2.get_arc(i).first, break2.get_arc(i).second);
    } 
  }

  if (break2.get_arc(0).first != Infty || break2.get_arc(1).first != Infty) {
    add_edge(index, break2.get_arc(0).first, break2.get_arc(1).first);
  }

  if (break2.get_arc(0).second != Infty || break2.get_arc(1).second != Infty) {
    add_edge(index, break2.get_arc(0).second, break2.get_arc(1).second);
  }
}

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_single_two_break(const TwoBreak<mcolor_t>& break2, partgraph_t& SG) { 
  for(size_t i = 0; i < 2; ++i) {
    if (break2.get_arc(i).first != Infty || break2.get_arc(i).second != Infty) {
      SG.erase(break2.get_arc(i).first, break2.get_arc(i).second);
    }
  }
	
  if (break2.get_arc(0).first != Infty || break2.get_arc(1).first != Infty) {
    SG.insert(break2.get_arc(0).first, break2.get_arc(1).first);
  }

  if (break2.get_arc(0).second != Infty || break2.get_arc(1).second != Infty) {
    SG.insert(break2.get_arc(0).second, break2.get_arc(1).second);
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_ins_del(const InsDel<mcolor_t>& insdel, bool record) {
  if (record) {
    insdel_history.push_back(insdel);
  }
 
  for (auto ic = insdel.cbegin_mcolor(); ic != insdel.cend_mcolor(); ++ic) {
	if (insdel.is_deletion_oper()) {
		erase_edge(ic->first, insdel.get_edge().first, insdel.get_edge().second);
	} else {
		add_edge(ic->first, insdel.get_edge().first, insdel.get_edge().second);
	} 
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_tandem_duplication(const TandemDuplication<mcolor_t>& dupl, bool record) {
  if (record) {
    tandem_dupl_history.push_back(dupl);
  }
 
  for (auto it = dupl.cbegin_edges(); it != dupl.cend_edges(); ++it) {
    for (auto im = dupl.cbegin_mcolor(); im != dupl.cend_mcolor(); ++im) {
	if (dupl.is_deletion_oper()) {
		erase_edge(im->first, it->first, it->second);
	} else {
		add_edge(im->first, it->first, it->second);
	}
    }
  }
} 

// check if a 2-break creates a circular chromosome
/*template<class graph_t, class mcolor_t>
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
}*/

#endif
