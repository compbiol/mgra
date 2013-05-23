/* 
** Module: Multicolors and Mutiple Breakpoint Graphs support
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

#ifndef MPBGRAPH_H_
#define MPBGRAPH_H_

#include <fstream>
#include <list>
#include <set>
#include <unordered_set>

#include "graph_colors.h"

#include "2break.h"
#include "mularcs.h"

#define member(S,x) ((S).find(x)!=(S).end())

extern std::ofstream outlog;

struct MBGraph {
  MBGraph(const std::vector<Genome>& genomes) { 
    local_graph.resize(genomes.size());
    std::unordered_set<orf_t> blocks; 

    for(size_t i = 0; i < genomes.size(); ++i) { 
      for(auto it = genomes[i].cbegin(); it != genomes[i].cend(); ++it) {
	if (blocks.count(it->second) == 0) { 
	  obverse_edges.insert(it->second + "t", it->second + "h");
	  blocks.insert(it->second);
	  vertex_set.insert(it->second + "t"); 
	  vertex_set.insert(it->second + "h"); 
	} 
      }
    } 
	
    for(size_t i = 0; i < genomes.size(); ++i) {
      add_edges(i, genomes[i], blocks);
    }	
  } 

  inline void add_edge(size_t index, const vertex_t& first, const vertex_t& second) { 
    local_graph[index].insert(first, second);
  }

  inline void erase_edge(size_t index, const vertex_t& first, const vertex_t& second) { 
    return local_graph[index].erase(first, second);
  } 

	inline void insert_twobreak(const TwoBreak<MBGraph, Mcolor>& break2) {
	 	history.push_back(break2);
	} 

	inline std::list<TwoBreak<MBGraph, Mcolor> > get_history() const { //FIXME: DEL
		return history; 
	} 

  //template<class mularcs_t, class colors_t>
  Mularcs<Mcolor> get_adjacent_multiedges(const vertex_t& u, const Graph_with_colors<Mcolor>& colors, const bool split_bad_colors = false) const; 
	
  inline size_t size_graph() const { 
    return local_graph.size();
  } 

	inline bool is_there_edge(size_t index, const vertex_t& first) const { 
		return local_graph[index].defined(first);		
	}

  inline vertex_t get_adj_vertex(const vertex_t& v) const {
    assert(obverse_edges.find(v) != obverse_edges.end());
    return obverse_edges.find(v)->second;
  } 
 
	inline vertex_t get_adj_vertex(size_t index, const vertex_t& first) const { 
		return local_graph[index][first];
	} 	 
	
        inline partgraph_t get_local_graph(size_t index) const { 
	  return local_graph[index];
	} 
	
	inline partgraph_t get_obverce_graph() const {
		return obverse_edges;
	}  
	
  inline std::set<vertex_t>::const_iterator begin_vertices() const { 
    return vertex_set.cbegin();
  } 
	
  inline std::set<vertex_t>::const_iterator end_vertices() const { 
    return vertex_set.cend();
  }

  inline std::vector<partgraph_t>::const_iterator begin_local_graphs() const { 
    return local_graph.cbegin();
  } 
	
  inline std::vector<partgraph_t>::const_iterator end_local_graphs() const { 
    return local_graph.cend(); 
  } 

	inline std::list<TwoBreak<MBGraph, Mcolor> >::const_iterator begin_history() const { 
		return history.cbegin();
	} 
	
	inline std::list<TwoBreak<MBGraph, Mcolor> >::const_iterator end_history() const { 
		return history.cend(); 
	} 

private:
	void add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks);
private:
	std::set<vertex_t> vertex_set;  // set of vertices //hash set? 
	partgraph_t obverse_edges; //OBverse relation 
	std::vector<partgraph_t> local_graph; // local graphs of each color 
	std::list<TwoBreak<MBGraph, Mcolor> > history;
};

typedef std::list<TwoBreak<MBGraph, Mcolor> > transform_t;
	
#endif
