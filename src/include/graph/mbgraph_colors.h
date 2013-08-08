/* 
** Module: Mutiple Breakpoint Graph with vec-T-consistent and T-consistent multicolors support
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

#ifndef MBGRAPH_COLOR_H_
#define MBGRAPH_COLOR_H_

#include "mbgraph.h"
#include "mularcs.h"

#include "genome_match.h" //FIXME REMOVE LATER

template<class mcolor_t>
struct mbgraph_with_colors: public MBGraph { 
  typedef typename std::set<mcolor_t>::const_iterator citer; 

  mbgraph_with_colors(const std::vector<Genome>& genomes, const ProblemInstance<Mcolor>& cfg); 

  bool is_simple_vertex(const vertex_t& v) const;
  bool is_indel_vertex(const vertex_t& v) const;  
  bool is_duplication_vertex(const vertex_t& v) const;
  bool is_have_self_loop(const vertex_t& v) const;

  Mularcs<mcolor_t> get_adjacent_multiedges(const vertex_t& u, bool split_bad_colors = false) const; 
  std::set<mcolor_t> split_color(const mcolor_t& color) const;
  bool are_adjacent_branches(const mcolor_t& A, const mcolor_t & B) const;

  inline mcolor_t get_complete_color() const {
    return complete_color;
  }

  inline mcolor_t get_complement_color(const mcolor_t& color) { 
    assert(color.is_one_to_one_match()); 
    if (CColorM.find(color) == CColorM.end()) { 
      mcolor_t temp = compute_complement_color(color);  
      CColorM.insert(std::make_pair(color, temp));
      CColorM.insert(std::make_pair(temp, color));
    }  
    return CColorM.find(color)->second;
  } 
  
  inline mcolor_t get_min_complement_color(const mcolor_t& color) {
    mcolor_t temp = get_complement_color(color);
    if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
      return color;
    } 
    return temp;
  }	

  inline bool is_T_consistent_color(const mcolor_t& color) const { 
    return (all_T_color.count(color) > 0);
  } 

  inline size_t count_vec_T_consitent_color() const { 
    return DiColor.size();
  } 
  
  inline bool is_vec_T_consistent_color(const mcolor_t& color) const {
    return (DiColor.count(color) > 0);
  }

  inline citer cbegin_T_consistent_color() const { 
    return DiColor.cbegin(); 
  } 

  inline citer cend_T_consistent_color() const { 
    return DiColor.cend(); 
  } 
private: 
  inline mcolor_t compute_complement_color(const mcolor_t& color) const {
      mcolor_t answer; 
      for(size_t j = 0; j < count_local_graphs(); ++j) { 
	if (!color.defined(j)) { 
	  answer.insert(j);
	} 
      } 
      return answer;
  }

private:
  mcolor_t complete_color;
  std::map<mcolor_t, mcolor_t> CColorM; //complementary multicolor
  std::set<mcolor_t> all_T_color; //all T-consistent colors
  std::set<mcolor_t> DiColor; // directed colors
}; 

template<class mcolor_t>
mbgraph_with_colors<mcolor_t>::mbgraph_with_colors(const std::vector<Genome>& genomes, const ProblemInstance<Mcolor>& cfg) 
: MBGraph(genomes)
{
  for (size_t i = 0; i < genomes.size(); ++i) {
    complete_color.insert(i);
    DiColor.insert(mcolor_t(i));
  }

  for (auto it = cfg.cbegin_trees(); it != cfg.cend_trees(); ++it) {
    it->build_vec_T_consistent_colors(DiColor);
  }

  DiColor.erase(complete_color);
  DiColor.erase(cfg.get_target());
  
  //check consistency
  for (auto id = DiColor.cbegin(); id != DiColor.cend(); ++id) {
    for(auto jd = id; jd != DiColor.end(); ++jd) {
      mcolor_t C(*id, *jd, mcolor_t::Intersection);
      if (!C.empty() && C.size() != id->size() && C.size() != jd->size()) {
	DiColor.erase(jd++);
	--jd;
      }
    }
  }
   
  for (const auto &vtc : DiColor) {
    all_T_color.insert(vtc);
    all_T_color.insert(compute_complement_color(vtc));
  }
}

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_simple_vertex(const vertex_t& v) const {
  Mularcs<mcolor_t> mularcs = get_adjacent_multiedges(v);
  if (mularcs.size() == 2 && mularcs.cbegin()->second.is_one_to_one_match() && mularcs.crbegin()->second.is_one_to_one_match() 
      && !is_duplication_vertex(v) && !is_indel_vertex(v)) {
	return true;  
  }
  return false; 
}  

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_have_self_loop(const vertex_t& v) const {
  Mularcs<mcolor_t> mularcs = get_adjacent_multiedges(v);
  if (mularcs.defined(v)) {
     return true;
  } 
  return false;
} 
 
template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_indel_vertex(const vertex_t& v) const {
  if (is_duplication_vertex(v)) {
    return false; 
  }  
 
  mcolor_t un = get_adjacent_multiedges(v).union_multicolors();
	
  if (!un.is_one_to_one_match() || (un == complete_color)) {
	return false; 
  }  

  return true; 
}

template<class mcolor_t>  
bool mbgraph_with_colors<mcolor_t>::is_duplication_vertex(const vertex_t& v) const {	
  if (is_have_self_loop(v)) {
    return true;
  } 

  Mularcs<mcolor_t> mularcs = get_adjacent_multiedges(v);
  for(auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) { 
    if (!im->second.is_one_to_one_match()) {
	return true; 
    }

    for(auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
      if (*im != *it) { 
        mcolor_t color(im->second, it->second, mcolor_t::Intersection);
        if (!color.empty()) { 
	  return true; 
        }
      }
    } 
  }  

  return false; 
} 


template<class mcolor_t>
Mularcs<mcolor_t> mbgraph_with_colors<mcolor_t>::get_adjacent_multiedges(const vertex_t& u, bool split_bad_colors) const { 
  if (u == Infty) {
    std::cerr << "mularcs ERROR: Infinite input" << std::endl;
    exit(1);
  }

  Mularcs<mcolor_t> output;
  for (size_t i = 0; i < count_local_graphs(); ++i) {
    if (local_graph[i].defined(u)) { 
      std::pair<partgraph_t::const_iterator, partgraph_t::const_iterator> iters = local_graph[i].equal_range(u);
      for (auto it = iters.first; it != iters.second; ++it) { 
	output.insert(it->second, i); 
      }
    } 
  }
  
  if (split_bad_colors) { 
    Mularcs<mcolor_t> split; 
    for(const auto &arc : output) {
      if (!is_vec_T_consistent_color(arc.second) && arc.second.size() < count_local_graphs()) {
	auto colors = split_color(arc.second);
	for(const auto &color : colors) {
	  split.insert(arc.first, color); 
	}
      } else { 
	split.insert(arc.first, arc.second); 
      }
    }
    return split; 
  }

  return output;
} 

/*
SplitColor(Q) представляет Q в виде дизъюнктного объединения T-consistent мультицветов, т.е. Q = Q1 U ... U Qm
где каждый Qi является T-consistent и все они попарно не пересекаются. SplitColor(Q) возвращает множество { Q1, Q2, ..., Qm }
(в частности, когда Q является T-consistent, имеем m=1 и Q1=Q).
Теперь, когда SplitBadColors = true, то и ребро (x,y) имеет мультицвет Q, то MBG.mulcols(x) будет содежать вместо (Q,x) пары:
(Q1,y), (Q2,y), ..., (Qm,y)
*/
template<class mcolor_t>
std::set<mcolor_t> mbgraph_with_colors<mcolor_t>::split_color(const mcolor_t& color) const {
  std::set<mcolor_t> S;

  if (is_vec_T_consistent_color(color)) {
    S.insert(color);
    return S;
  }

  equivalence<size_t> EQ;
  for(auto iq = color.cbegin(); iq != color.cend(); ++iq) { 
    EQ.addrel(iq->first, iq->first);
  } 

  for (const auto &vtc: DiColor) { 
    mcolor_t C(vtc, color, mcolor_t::Intersection);
    if (C.size() >= 2 && C.size() == vtc.size() ) {
      for (auto iq = C.cbegin(); iq != C.cend(); ++iq) {
	EQ.addrel(iq->first, C.cbegin()->first);
      }
    }
  }

  EQ.update();
  std::map<size_t, mcolor_t> cls = EQ.get_eclasses<mcolor_t>(); 
  for(const auto &col : cls) {
    S.insert(col.second);
  }
  return S;
} 

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::are_adjacent_branches(const mcolor_t& A, const mcolor_t & B) const {
  if (!is_T_consistent_color(A) || !is_T_consistent_color(B)) { 
    return false;
  } 

  mcolor_t Q1; 
  mcolor_t Q2;
  
  if (A.size() >= B.size()) {
    Q1 = A;
    Q2 = B;
  } else {
    Q1 = B;
    Q2 = A;
  }
  
  mcolor_t C(Q1, Q2, mcolor_t::Difference);
  
  if (C.size() == Q1.size() - Q2.size() && is_T_consistent_color(C)) { 		
    return true;
  } 
  
  mcolor_t M(Q1, Q2, mcolor_t::Union); 
  if (M.size() == Q1.size() + Q2.size() && is_T_consistent_color(M)) { 	
    return true;
  } 
  
  return false;
}

#endif

