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

#include "defined.h" 
#include "genome.h"

struct MBGraph {
  typedef std::string orf_t;

  MBGraph(const std::vector<Genome>& genomes) 
  : local_graph(genomes.size()) 
  { 
    std::unordered_set<orf_t> blocks;

    for (const auto& genome : genomes) { 
      for(const auto &chromosome : genome) {
	for(const auto &orf : chromosome.second) {
	  if (blocks.count(orf.second.first) == 0) { 
	    obverse_edges.insert(orf.second.first + "t", orf.second.first + "h");
	    blocks.insert(orf.second.first);
	    vertex_set.insert(orf.second.first + "t"); 
	    vertex_set.insert(orf.second.first + "h"); 
	  }
        } 
      }
    }

    for(size_t i = 0; i < genomes.size(); ++i) {
      add_edges(i, genomes[i]); 
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

  inline vertex_t get_obverse_vertex(const vertex_t& v) const {
    assert(obverse_edges.count(v) != 0);
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
    assert(index < local_graph.size() && (local_graph[index].count(first) != 0));
    return local_graph[index].find(first)->second;
  } 	 

  inline size_t size() const { 
    return vertex_set.size(); 
  }
 
  inline size_t count_local_graphs() const { 
    return local_graph.size();
  } 

  inline partgraph_t get_local_graph(size_t index) const {
    return local_graph[index];
  } 

  inline obverse_graph_t get_obverse_graph() const {
     return obverse_edges;
  }
 
  inline std::set<vertex_t>::iterator begin() { 
    return vertex_set.begin();
  } 
	
  inline std::set<vertex_t>::iterator end() { 
    return vertex_set.end();
  }

  inline std::set<vertex_t>::const_iterator begin() const { 
    return vertex_set.cbegin();
  } 
	
  inline std::set<vertex_t>::const_iterator end() const { 
    return vertex_set.cend();
  }
		
  inline std::vector<partgraph_t>::const_iterator cbegin_local_graphs() const { 
    return local_graph.cbegin();
  } 
	
  inline std::vector<partgraph_t>::const_iterator cend_local_graphs() const { 
    return local_graph.cend(); 
  } 

private:
  void add_edges(size_t index, const Genome& genome) {
    auto rearLambda = [] (const std::pair<orf_t, int> & orf) -> vertex_t { 
      return ((orf.second > 0)?(orf.first + "h"):(orf.first + "t"));
    };

    auto frontLambda = [] (const std::pair<orf_t, int> & orf) -> vertex_t { 
      return ((orf.second > 0)?(orf.first + "t"):(orf.first + "h"));
    };

    for (const auto &chromosome : genome) {
      vertex_t current_vertex = rearLambda(chromosome.second.begin()->second);
      for (auto gene = (++chromosome.second.begin()); gene != chromosome.second.end(); ++gene) { 
        local_graph[index].insert(current_vertex, frontLambda(gene->second));	
        current_vertex = rearLambda(gene->second);
      }

      if (chromosome.second.is_circular()) {
	local_graph[index].insert(frontLambda(chromosome.second.begin()->second), rearLambda((--chromosome.second.end())->second)); 
      } else { 
	local_graph[index].insert(frontLambda(chromosome.second.begin()->second), Infty);
	local_graph[index].insert(rearLambda((--chromosome.second.end())->second), Infty); 
      }
    }
  }

protected:		
  std::set<vertex_t> vertex_set; //set of vertice
  obverse_graph_t obverse_edges; //obverse relation 
  std::vector<partgraph_t> local_graph; //local graphs of each color 
};	

#endif
