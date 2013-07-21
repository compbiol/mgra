/* 
** Module: Mutiple Breakpoint Graph with vec-T-consistent and T-consistent multicolors and save all events in evolution support
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008 - 2013 by Max Alekseyev <maxal@cse.sc.edu> and Pavel Avdeyev
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

#define member(S,x) ((S).find(x)!=(S).end())

template<class mcolor_t>
struct mbgraph_with_history : public mbgraph_with_colors<mcolor_t> { 

	mbgraph_with_history(const std::vector<Genome>& genomes, const ProblemInstance& cfg)
	: mbgraph_with_colors<mcolor_t>(genomes, cfg) 
	{ 
	} 

	void apply_two_break(const TwoBreak<mcolor_t>& break2, bool record = false);
	void apply_single_two_break(size_t index, const TwoBreak<mcolor_t>& break2);  
	void apply_single_two_break(const TwoBreak<mcolor_t>& break2, partgraph_t& SG); //FIXME DEL 

	inline typename std::list<TwoBreak<mcolor_t> > get_history() const { //FIXME: DEL
		return history; 
	} 

	inline typename std::list<TwoBreak<mcolor_t> >::const_iterator begin_history() const { 
		return history.cbegin();
	} 
	
	inline typename std::list<TwoBreak<mcolor_t> >::const_iterator end_history() const { 
		return history.cend(); 
	} 
private: 
	std::list<TwoBreak<mcolor_t> > history;	
};

typedef std::list<TwoBreak<Mcolor> > transform_t;

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_two_break(const TwoBreak<mcolor_t>& break2, bool record) { 
  if (record) {
    history.push_back(break2);
  }
 
  for (auto ic = break2.MultiColor.cbegin(); ic != break2.MultiColor.cend(); ++ic) {
    for (size_t i = 0; i < 2; ++i) {
      if (break2.OldArc[i].first != Infty && break2.OldArc[i].second != Infty) {
	erase_edge(ic->first, break2.OldArc[i].first, break2.OldArc[i].second);
      } 
    }

    if (break2.OldArc[0].first != Infty && break2.OldArc[1].first != Infty) {
      add_edge(ic->first, break2.OldArc[0].first, break2.OldArc[1].first);
    }

    if (break2.OldArc[0].second != Infty && break2.OldArc[1].second != Infty) {
      add_edge(ic->first, break2.OldArc[0].second, break2.OldArc[1].second);
    }
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_single_two_break(size_t index, const TwoBreak<mcolor_t>& break2) {
  for(size_t i = 0; i < 2; ++i) {
    if (break2.OldArc[i].first != Infty && break2.OldArc[i].second != Infty) {
       erase_edge(index, break2.OldArc[i].first, break2.OldArc[i].second);
    }
  }

  if (break2.OldArc[0].first != Infty && break2.OldArc[1].first != Infty) {
    add_edge(index, break2.OldArc[0].first, break2.OldArc[1].first);
  }

  if (break2.OldArc[0].second != Infty && break2.OldArc[1].second != Infty) {
    add_edge(index, break2.OldArc[0].second, break2.OldArc[1].second);
  }
}

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_single_two_break(const TwoBreak<mcolor_t>& break2, partgraph_t& SG) { 
  for(size_t i = 0; i < 2; ++i) {
    if (break2.OldArc[i].first != Infty && break2.OldArc[i].second != Infty) {
      SG.erase(break2.OldArc[i].first, break2.OldArc[i].second);
    }
  }
	
  if (break2.OldArc[0].first != Infty && break2.OldArc[1].first != Infty) {
    SG.insert(break2.OldArc[0].first, break2.OldArc[1].first);
  }

  if (break2.OldArc[0].second != Infty && break2.OldArc[1].second != Infty) {
    SG.insert(break2.OldArc[0].second, break2.OldArc[1].second);
  }
} 

#endif
