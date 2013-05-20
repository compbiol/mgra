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
//#include "utility/sym_multi_hashmap.h"

#include "2break.h"
#include "mularcs.h"


#define member(S,x) ((S).find(x)!=(S).end())

extern std::ofstream outlog;
//const vertex_t Infty = "oo";

typedef sym_multi_hashmap<vertex_t> partgraph_t;

struct MBGraph {

	inline MBGraph(const std::vector<Genome>& genomes) { 
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

	inline void insert_twobreak(const TwoBreak<MBGraph>& break2) {
	 	history.push_back(break2);
	} 

	inline std::list<TwoBreak<MBGraph> > get_history() const { //FIXME: DEL
		return history; 
	} 

	//template<class mularcs_t, class colors_t>
	Mularcs get_adjacent_multiedges(const vertex_t& u, const ColorsGraph<Mcolor>& colors, const bool split_bad_colors = false) const; 
	
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

	inline std::list<TwoBreak<MBGraph> >::const_iterator begin_history() const { 
		return history.cbegin();
	} 
	
	inline std::list<TwoBreak<MBGraph> >::const_iterator end_history() const { 
		return history.cend(); 
	} 

private:
	std::set<Mcolor> split_color(const Mcolor& Q, const ColorsGraph<Mcolor>& colors, bool split_bad_color) const; //FIXME: go to Mcolor	
	void add_edges(size_t index, const Genome& genome, const std::unordered_set<orf_t>& blocks);

private:
	std::set<vertex_t> vertex_set;  // set of vertices //hash set? 
	partgraph_t obverse_edges; //OBverse relation 
	std::vector<partgraph_t> local_graph; // local graphs of each color 
	std::list<TwoBreak<MBGraph> > history;
};

typedef std::list<TwoBreak<MBGraph> > transform_t;
//extern transform_t History;



/*template<class mularcs_t, class colors_t>
mularcs_t MBGraph::get_adjacent_multiedges(const vertex_t& u, const colors_t& colors, const bool split_bad_colors) const {
	if (u == Infty) {
		std::cerr << "mularcs ERROR: Infinite input" << std::endl;
		exit(1);
	}

	mularcs_t output;
	for (int i = 0; i < size_graph(); ++i) {
		if (local_graph[i].defined(u)) { 
			std::pair<partgraph_t::const_iterator, partgraph_t::const_iterator> iters = local_graph[i].equal_range(u);
			for (auto it = iters.first; it != iters.second; ++it) { 
				if (output.find(it->second) != output.cend()) { 
					output.find(it->second)->second.insert(i);
				} else { 
					output.insert(it->second, Mcolor(i));	
				} 
			}
		} else { 
			if (output.find(Infty) != output.cend()) { 
				output.find(Infty)->second.insert(i);
			} else { 
				output.insert(Infty, Mcolor(i));
			} 
		} 
	}

	if (split_bad_colors) { 
		mularcs_t split; 
		for(auto im = output.cbegin(); im != output.cend(); ++im) {
			if (!member(colors.DiColor, im->second) && im->second.size() < size_graph()) {
				auto C = split_color(im->second, colors, split_bad_colors);
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
} */
	
#endif
