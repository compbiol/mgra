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
#include "utility/sym_map.h"

#include "genome_match.h" //FIXME REMOVE LATER

template<class mcolor_t>
struct mbgraph_with_colors: public MBGraph { 
  typedef typename std::set<mcolor_t>::const_iterator citer; 

  mbgraph_with_colors(const std::vector<Genome>& genomes, const ProblemInstance<Mcolor>& cfg); 

  void update_complement_color(const std::vector<mcolor_t>& colors); //FIXME: THINK ABOUT IT

  bool is_simple_vertex(const vertex_t& v) const;
  bool is_indel_vertex(const vertex_t& v) const;  
  bool is_duplication_vertex(const vertex_t& v) const;
  bool is_have_self_loop(const vertex_t& v) const;

  Mularcs<mcolor_t> get_adjacent_multiedges(const vertex_t& u, bool split_bad_colors = false) const; 
  bool are_adjacent_branches(const mcolor_t& A, const mcolor_t & B) const;

  inline mcolor_t get_complete_color() const {
    return complete_color;
  }

  inline mcolor_t get_complement_color(const mcolor_t& color) const { 
    assert (CColorM.find(color) != CColorM.end());
    return CColorM.find(color)->second;
  } 
  
  inline mcolor_t get_min_complement_color(const mcolor_t& color) const {
    mcolor_t temp = get_complement_color(color);
    if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
      return color;
    } 
    return temp;
  }	

  inline bool is_T_consistent_color(const mcolor_t& color) const { 
    return (all_T_color.count(color) > 0);
  } 

  inline size_t count_vec_T_color() const { 
    return DiColor.size();
  } 
  
  inline bool is_vec_T_color(const mcolor_t& color) const {
	return (DiColor.count(color) > 0);
  }

  inline citer cbegin_T_color() const { 
	return DiColor.cbegin(); 
  } 

  inline citer cend_T_color() const { 
	return DiColor.cend(); 
  } 
private: 
  inline mcolor_t compute_complement_color(const mcolor_t& color) const {
      mcolor_t answer; 
      for(size_t j = 0; j < size_graph(); ++j) { 
	if (!color.mymember(j)) { 
	  answer.insert(j);
	} 
      } 
      return answer;
  }
private:
  mcolor_t complete_color;
  sym_map<mcolor_t> CColorM; //complementary multicolor
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
    it->get_dicolors(DiColor);
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
   
  for (auto id = DiColor.begin(); id != DiColor.end(); ++id) {
    all_T_color.insert(*id);
    all_T_color.insert(compute_complement_color(*id));
  }
}

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_simple_vertex(const vertex_t& v) const {
  Mularcs<mcolor_t> mularcs = get_adjacent_multiedges(v);
  if (mularcs.size() == 2) { 
    if (mularcs.cbegin()->second.is_one_to_one_match() && mularcs.crbegin()->second.is_one_to_one_match()) { 
      if (!is_duplication_vertex(v) && !is_indel_vertex(v)) {
	return true;  
      }
    } 
  } 
  return false; 
}  

template<class mcolor_t>
bool mbgraph_with_colors<mcolor_t>::is_have_self_loop(const vertex_t& v) const {
  Mularcs<mcolor_t> mularcs = get_adjacent_multiedges(v);
  if (mularcs.find(v) != mularcs.cend()) {
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
	
  if (!un.is_one_to_one_match()) {
	return false; 
  }  

  if (un == complete_color) { 
    return false;
  }
 
  return true; 
}

template<class mcolor_t>  
bool mbgraph_with_colors<mcolor_t>::is_duplication_vertex(const vertex_t& v) const {	
  if (this->is_have_self_loop(v)) {
    return true;
 } 

  Mularcs<mcolor_t> mularcs = get_adjacent_multiedges(v);
  for(auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) { 
    for(auto it = mularcs.cbegin(); it != mularcs.cend(); ++it) {
      if (*im == *it) { 
	continue;
      } 

      mcolor_t color(im->second, it->second, mcolor_t::Intersection);
      if (!color.empty()) { 
	return true; 
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
  for (size_t i = 0; i < size_graph(); ++i) {
    if (local_graph[i].defined(u)) { 
      std::pair<partgraph_t::const_iterator, partgraph_t::const_iterator> iters = local_graph[i].equal_range(u);
      for (auto it = iters.first; it != iters.second; ++it) { 
	if (output.find(it->second) != output.cend()) { 
	  output.find(it->second)->second.insert(i);
	} else { 
	  output.insert(it->second, mcolor_t(i));	
	} 
      }
    } 
  }
  
  if (split_bad_colors) { 
    Mularcs<mcolor_t> split; 
    for(auto im = output.cbegin(); im != output.cend(); ++im) {
      if (!is_vec_T_color(im->second) && im->second.size() < size_graph()) {
	auto C = im->second.split_color(*this, split_bad_colors);
	for(auto ic = C.begin(); ic != C.end(); ++ic) {
	  split.insert(im->first, *ic); 
	}
      } else { 
	split.insert(im->first, im->second); 
      }
    }
    return split; 
  }

  return output;
} 

template<class mcolor_t>
void mbgraph_with_colors<mcolor_t>::update_complement_color(const std::vector<mcolor_t>& colors) {
  for(auto it = colors.begin(); it != colors.end(); ++it) { 
    if (!CColorM.defined(*it)) { 
      mcolor_t temp = compute_complement_color(*it);  
      CColorM.insert(*it, temp);
    } 
  } 
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

