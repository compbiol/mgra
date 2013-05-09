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

#include "genome_match.h"
#include "structures/Tree.h"

#include "utility/equivalence.h"
#include "utility/sym_multi_hashmap.h"
#include "utility/sym_hashmap.h"
#include "utility/sym_map.h"

#define member(S,x) ((S).find(x)!=(S).end())

extern std::ofstream outlog;
const std::string Infty = "oo";

typedef std::string vertex_t;

typedef sym_multi_hashmap<vertex_t> multi_hashmap;

typedef sym_multi_hashmap<vertex_t> partgraph_t; //s graph represented in list of edges/ but not

typedef std::map <vertex_t, Mcolor> mularcs_t;
typedef std::multimap <vertex_t, Mcolor> multimularcs_t;

/*Unordered set of vertex is bad, because need order. Why? */
struct MBGraph {
	void init(const std::vector<Genome>& genome, const ProblemInstance& cfg); //it's constructor

	/*function for colors*/	
	void update_complement_color(const std::vector<Mcolor>& colors);

	/*inline bool is_have_complement_color(const Mcolor& color) const { 
		CColorM.find(color) != CColorM.end()
	} */

	inline Mcolor get_complement_color(const Mcolor& color) const { 
		assert(CColorM.find(color) != CColorM.end());
		return CColorM.find(color)->second;
	} 

	inline Mcolor get_min_complement_color(const Mcolor& color) const {
		Mcolor temp = get_complement_color(color);
		if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
			return color;
		} 
		return temp;
	}	
	
	inline bool is_T_consistent_color(const Mcolor& col) const { 
		return (all_T_color.find(col) != all_T_color.end());
	} 

	/*function for graphs*/
	mularcs_t get_adjacent_multiedges(const vertex_t& u) const; 
	//multimularcs_t get_adjacent_multiedges_v2(const vertex_t& u) const;
	multimularcs_t get_adjacent_multiedges_with_split(const vertex_t& u) const; 

	std::map<std::pair<Mcolor, Mcolor>, size_t> get_count_Hsubgraph() const; 

	inline bool is_simple_vertice(const mularcs_t& adj_edges) const {
		if (adj_edges.size() == 2 && adj_edges.begin()->second.is_good_multiedge() && adj_edges.rbegin()->second.is_good_multiedge()) { 
			return true; 
		} 
		return false; 
	}  

	inline bool is_fair_vertice(const mularcs_t& adj_edges) const {
		if (adj_edges.size() == 3) {
			for(auto it = adj_edges.cbegin(); it != adj_edges.cend(); ++it) {
				if (!it->second.is_good_multiedge()) { 
					return false;
				}  
			} 
			return true; 
		}  
		return false; 
	} 

	inline size_t size_graph() const { 
		return LG.size();
	} 

	inline partgraph_t get_obverce_graph() const { //it's big function. and use only gen_dist 
		return obverse_edges;
	}  

	inline vertex_t get_adj_vertex(const vertex_t& v) const {
		assert(obverse_edges.find(v) != obverse_edges.end());
		return obverse_edges.find(v)->second;
	} 

 	inline std::set<std::string>::const_iterator begin_vertices() const { 
		return vertex_set.cbegin();
	} 
	
	inline std::set<std::string>::const_iterator end_vertices() const { 
		return vertex_set.cend();
	} 

	inline partgraph_t::const_iterator begin_partgraph(size_t index) const { 
		return LG[index].cbegin();
	} 
	
	inline partgraph_t::const_iterator end_partgraph(size_t index) const { 
		return LG[index].cend(); 
	} 

	/************Other**************/
	inline void erase_vertex_in_graph(const vertex_t& vertex) { 
		obverse_edges.erase(vertex);
	} 
private:
	void parsing_tree(const std::vector<Genome>& genomes, const ProblemInstance& cfg);
	Mcolor add_tree(const std::string& tree, std::vector<std::string>& output);

	std::set<Mcolor> split_color(const Mcolor& Q) const;
	
	partgraph_t add_edges(const Genome& genome, const std::unordered_set<orf_t>& blocks) const;
	void build_graph(const std::vector<Genome>& genomes); 	
private:
	//edges
	partgraph_t obverse_edges; //OBverse relation 
	std::set<std::string> vertex_set;  // set of vertices //rename and take private and hash set

	//color
	sym_map<Mcolor> CColorM; //complementary multicolor
	std::set<Mcolor> all_T_color;	//all T-consistent colors	
public:
	std::vector<partgraph_t> LG; // local graphs of each color //rename and take private 	
	
	std::set<Mcolor> DiColor;   // directed colors
	std::vector<Mcolor> TColor; // colors corresponding to ancestral genomes. nodes in tree
	std::set<Mcolor> DiColorUnsafe; //stage 4 wtf

	static bool SplitBadColors;

	const Mcolor& CColor(const Mcolor& S);

	// complementary colors representative
	Mcolor CColorRep(const Mcolor& c);

	bool AreAdjacentBranches(const Mcolor&, const Mcolor&) const;
};

extern MBGraph MBG; //deel
#endif
