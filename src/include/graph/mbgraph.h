/* 
** Module: Mutiple Breakpoint Graph support
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

#ifndef MBGRAPH_H_
#define MBGRAPH_H_

#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <cassert>

#include "utility/sym_multi_hashmap.h"
#include "2break.h"
#include "genome.h"
#include "2break.h"
//typedef std::string vertex_t;
//typedef sym_multi_hashmap<vertex_t> partgraph_t;

//const vertex_t Infty = "oo"; 

struct MBGraph {
  MBGraph(const std::vector<Genome>& genomes) 
  : local_graph(genomes.size()) 
  { 
    std::unordered_set<orf_t> blocks; 

    for(size_t i = 0; i < genomes.size(); ++i) { 
      for(auto it = genomes[i].cbegin(); it != genomes[i].cend(); ++it) {
	if (blocks.count(it->second.first) == 0) { 
	  obverse_edges.insert(it->second.first + "t", it->second.first + "h");
	  blocks.insert(it->second.first);
	  vertex_set.insert(it->second.first + "t"); 
	  vertex_set.insert(it->second.first + "h"); 
	} 
      }
    } 

    for(size_t i = 0; i < genomes.size(); ++i) {
      add_edges(i, genomes[i], blocks);
    }	
  } 

  inline void add_edge(size_t index, const vertex_t& first, const vertex_t& second) { 
    assert(index < local_graph.size());
    local_graph[index].insert(first, second);
  }

  inline void erase_edge(size_t index, const vertex_t& first, const vertex_t& second) {
    assert(index < local_graph.size());
    return local_graph[index].erase(first, second);
  } 

  inline size_t size_graph() const { 
    return local_graph.size();
  } 

  inline size_t count_vertex() const { 
    return vertex_set.size();
  } 
	
  inline vertex_t get_obverse_vertex(const vertex_t& v) const {
    assert(obverse_edges.find(v) != obverse_edges.cend());
    return obverse_edges.find(v)->second;
  } 

  // FIXME: PROBLEM WITH TARGET. I NOT UNDERSTAND ALGORITHM ABOUT THIS. If changed, chaged this function
  inline bool is_exist_edge(size_t index, const vertex_t& first) const { 
    assert(index < local_graph.size());
    if (local_graph[index].find(first) != local_graph[index].end()) { 
	if (local_graph[index].find(first)->second == Infty)  {
		return false;
	}
	return true;
    } else { 
	return false; 
    }  
    //return local_graph[index].defined(first);		
  }
 
  //FIXME IF WE RECONSTRUCT ANCESTORS WITH DUPLICATION EQUAL RANGE
  inline vertex_t get_adjecent_vertex(size_t index, const vertex_t& first) const {  
    assert(index < local_graph.size());
    assert(local_graph[index].find(first) != local_graph[index].cend());
    return local_graph[index].find(first)->second;
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

private:
  void add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks) {
    vertex_t first_vertex; //in chromosome
    vertex_t current_vertex; // is rightmost vertex of the previous block
    vertex_t prev_chr;

    for(auto iter = genome.cbegin(); iter != genome.cend(); ++iter) {
      if (blocks.find(iter->second.first) != blocks.end()) { 
	if (iter->first.first == prev_chr) {
          if (iter->second.second > 0) {
	    local_graph[index].insert(current_vertex, iter->second.first + "t");
	  } else {
	    local_graph[index].insert(current_vertex, iter->second.first + "h");
	  }
	} else { // new chromosome detected				 
	  if (!first_vertex.empty()) { 
	    if (genome.isCircular(prev_chr)) {
		local_graph[index].insert(first_vertex, current_vertex);
	    } else { 
		local_graph[index].insert(first_vertex, Infty);
		local_graph[index].insert(current_vertex, Infty);
	    } 
	  }

	  if (iter->second.second > 0) { 
	    first_vertex = iter->second.first + "t"; 
	  } else { 
	    first_vertex = iter->second.first + "h";
	  } 

	  prev_chr = iter->first.first;
	} 
		
	if (iter->second.second > 0) { 
	  current_vertex = iter->second.first + "h"; 
	} else { 
	  current_vertex = iter->second.first + "t";
	} 
      } 
    }

    if (!first_vertex.empty()) { 
      if (genome.isCircular(prev_chr)) {
        local_graph[index].insert(first_vertex, current_vertex);
      } else { 
	local_graph[index].insert(first_vertex, Infty);
	local_graph[index].insert(current_vertex, Infty);
      }  
    } 
  }
protected:
  std::set<vertex_t> vertex_set;  // set of vertices //hash set? 
  partgraph_t obverse_edges; //OBverse relation 
  std::vector<partgraph_t> local_graph; // local graphs of each color 
};	

#endif
