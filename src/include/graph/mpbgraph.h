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
#include <sstream>
#include <list>
#include <map>
#include <set>
#include <unordered_set>
using namespace std;

#include "graph_colors.h"
#include "structures/Tree.h"

#include "utility/equivalence.h"
#include "utility/sym_multi_hashmap.h"

#include "mularcs.h"
 
#define member(S,x) ((S).find(x)!=(S).end())

extern std::ofstream outlog;
const std::string Infty = "oo";

typedef sym_multi_hashmap<vertex_t> partgraph_t;

//template<class mularcs_t>
struct MBGraph {
	MBGraph(const std::vector<Genome>& genome, const ProblemInstance& cfg);

	inline void add_edge(size_t index, const std::string& first, const std::string& second) { 
		local_graph[index].insert(first, second);
	}

	inline void erase_edge(size_t index, const std::string& first, const std::string& second) { 
		return local_graph[index].erase(first, second);
	} 

	Mularcs get_adjacent_multiedges(const vertex_t& u, const ColorsGraph<Mcolor>& colors, const bool split_bad_colors = false) const; 
	
	inline std::set<std::string>::const_iterator begin_vertices() const { 
		return vertex_set.cbegin();
	} 
	
	inline std::set<std::string>::const_iterator end_vertices() const { 
		return vertex_set.cend();
	}

	inline size_t size_graph() const { 
		return local_graph.size();
	} 

	inline bool is_there_edge(size_t index, const std::string& first) const { 
		return local_graph[index].defined(first);		
	}

	inline vertex_t get_adj_vertex(const vertex_t& v) const {
		assert(obverse_edges.find(v) != obverse_edges.end());
		return obverse_edges.find(v)->second;
	} 
 
	inline std::string get_adj_vertex(size_t index, const std::string& first) const { 
		return local_graph[index][first];
	} 	 
	
	inline partgraph_t get_local_graph(size_t index) const { 
		return local_graph[index];
	} 
	
	inline partgraph_t get_obverce_graph() const {
		return obverse_edges;
	}  
	
	inline std::vector<partgraph_t>::const_iterator begin_local_graphs() const { 
		return local_graph.cbegin();
	} 
	
	inline std::vector<partgraph_t>::const_iterator end_local_graphs() const { 
		return local_graph.cend(); 
	} 

private:
	std::set<Mcolor> split_color(const Mcolor& Q, const ColorsGraph<Mcolor>& colors, bool split_bad_color) const; //FIXME: go to Mcolor	
	void add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks);

private:
	std::set<vertex_t> vertex_set;  // set of vertices //hash set? 
	partgraph_t obverse_edges; //OBverse relation 
	std::vector<partgraph_t> local_graph; // local graphs of each color 
};

#endif
