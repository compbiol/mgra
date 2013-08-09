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

template<class mcolor_t>
struct mbgraph_with_colors: public MBGraph { 
  typedef typename std::set<mcolor_t>::const_iterator citer; 

  template<class pconf_t>
  mbgraph_with_colors(const std::vector<Genome>& genomes, const pconf_t& cfg); 

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
    if (compliment_colors.count(color) == 0) { 
      mcolor_t temp = compute_complement_color(color);  
      compliment_colors.insert(std::make_pair(color, temp));
      compliment_colors.insert(std::make_pair(temp, color));
      return temp;
    }  
    return compliment_colors.find(color)->second;
  } 
  
  inline mcolor_t get_min_complement_color(const mcolor_t& color) {
    mcolor_t temp = get_complement_color(color);
    if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
      return color;
    } 
    return temp;
  }	

  inline bool is_T_consistent_color(const mcolor_t& color) const { 
    return (T_consistent_colors.count(color) > 0);
  } 

  inline size_t count_vec_T_consitent_color() const { 
    return vec_T_consistent_colors.size();
  } 
  
  inline bool is_vec_T_consistent_color(const mcolor_t& color) const {
    return (vec_T_consistent_colors.count(color) > 0);
  }

  inline citer cbegin_T_consistent_color() const { 
    return vec_T_consistent_colors.cbegin(); 
  } 

  inline citer cend_T_consistent_color() const { 
    return vec_T_consistent_colors.cend(); 
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

protected:
  mcolor_t complete_color;
  std::map<mcolor_t, mcolor_t> compliment_colors;
  std::set<mcolor_t> T_consistent_colors;
  std::set<mcolor_t> vec_T_consistent_colors;
}; 

template<class mcolor_t>
template<class pconf_t>
mbgraph_with_colors<mcolor_t>::mbgraph_with_colors(const std::vector<Genome>& genomes, const pconf_t& cfg) 
: MBGraph(genomes)
{
  for (size_t i = 0; i < genomes.size(); ++i) {
    complete_color.insert(i);
    vec_T_consistent_colors.insert(mcolor_t(i));
  }

  for (auto it = cfg.cbegin_trees(); it != cfg.cend_trees(); ++it) {
    it->build_vec_T_consistent_colors(vec_T_consistent_colors);
  }

  vec_T_consistent_colors.erase(complete_color);
  vec_T_consistent_colors.erase(cfg.get_target());
  
  //check consistency
  for (auto id = vec_T_consistent_colors.cbegin(); id != vec_T_consistent_colors.cend(); ++id) {
    for(auto jd = id; jd != vec_T_consistent_colors.end(); ++jd) {
      mcolor_t color(*id, *jd, mcolor_t::Intersection);
      if (!color.empty() && color.size() != id->size() && color.size() != jd->size()) {
	vec_T_consistent_colors.erase(jd++);
	--jd;
      }
    }
  }
   
  for (const auto &vtc : vec_T_consistent_colors) {
    T_consistent_colors.insert(vtc);
    T_consistent_colors.insert(compute_complement_color(vtc));
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
  std::set<mcolor_t> answer;

  if (is_vec_T_consistent_color(color)) {
    answer.insert(color);
  } else { 
    equivalence<size_t> equiv;
    for(auto iq = color.cbegin(); iq != color.cend(); ++iq) { 
      equiv.addrel(iq->first, iq->first);
    } 

    for (const auto &vtc: vec_T_consistent_colors) { 
      mcolor_t inter_color(vtc, color, mcolor_t::Intersection);
      if (inter_color.size() >= 2 && inter_color.size() == vtc.size() ) {
        for (const auto &col : inter_color) { 
	 equiv.addrel(col.first, inter_color.cbegin()->first);
        }
      }
    }

    equiv.update();
    std::map<size_t, mcolor_t> classes = equiv.get_eclasses<mcolor_t>(); 
    for(const auto &col : classes) {
      answer.insert(col.second);
    }
  }
  return answer;
} 

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::are_adjacent_branches(const mcolor_t& mcolor_a, const mcolor_t & mcolor_b) const {
  if (!is_T_consistent_color(mcolor_a) || !is_T_consistent_color(mcolor_b)) { 
    return false;
  } 

  mcolor_t color1; 
  mcolor_t color2;
  
  if (mcolor_a.size() >= mcolor_b.size()) {
    color1 = mcolor_a;
    color2 = mcolor_b;
  } else {
    color1 = mcolor_b;
    color2 = mcolor_a;
  }
  
  mcolor_t diff_color(color1, color2, mcolor_t::Difference);
  
  if (diff_color.size() == color1.size() - color2.size() && is_T_consistent_color(diff_color)) { 		
    return true;
  } 
  
  mcolor_t union_color(color1, color2, mcolor_t::Union); 
  if (union_color.size() == color1.size() + color2.size() && is_T_consistent_color(union_color)) { 	
    return true;
  } 
  
  return false;
}

#endif

