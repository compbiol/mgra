#ifndef ESTIMATE_H_
#define ESTIMATE_H_

#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>

#include "mcolor.h"
#include "mpbgraph.h"

typedef std::string vertex_t;

struct Statistics { 
  Statistics(const MBGraph& gr); 

  std::vector<std::string> get_compl_stat() const;   
  std::vector<std::string> get_no_compl_stat() const;	
  std::vector<Mcolor> get_new_color() const;

  std::map<std::pair<Mcolor, Mcolor>, size_t> get_Hsubgraph(); //count H-subgraph for stage2 

private:
  void count_weak_simple_vertex(); 
  void count_compl_multiedges(); //count good edges for stage1
  void count_not_compl_multiedges(); 
  void count_cycles();
  void count_chromosomes();

  __attribute__((always_inline)) inline size_t calc_value(const std::map<Mcolor, size_t>& where, const Mcolor& what) const { 
    if (where.find(what) != where.end()) { 
      return where.find(what)->second;
    } 
    return 0; 
  } 

private: 
  const MBGraph& graph;

  //vertices
  std::unordered_map<size_t, size_t > multidegree_count; // multidegree_count[n] = # vertices of multidegree n. 
  std::map<Mcolor, size_t> simple_vertices_count;  	 // simple_vertices_count[min(S,!S)] = # simple vertices incident to S-colored
  std::map<Mcolor, size_t> simple_vertices_alone_count;  // simple_vertices_alone_count[min(S,!S)] = # simple vertices incident to S-colored, with no good neighbors

  //edges
  std::map<Mcolor, size_t> not_compl_multiedges_count;	// not complement multiedges[S] = # multiedges of not complement multicolor S,
  std::map<Mcolor, size_t> compl_multiedges_count; 	// multiedges_count[S] = # multiedges of multicolor S.
  std::map<Mcolor, size_t> good_multiedges_count; 	// good_multiedges_count[S] = # good multiedges of multicolor S. 
  std::map<Mcolor, size_t> good_irrer_multiedges_count;	// ME[S] = # good irregular multiedges of multicolor S.
  std::map<Mcolor, size_t> simple_multiedges_count;	// ME[S] = # simple multiedges of multicolor S.

  std::map<std::pair<Mcolor, Mcolor>, size_t> Hcount; // count H-subgraphs
  //std::map<std::pair<Mcolor, Mcolor>, bool> Hmid;   // middle edge is T-consistent?
	
  //cycles
  std::map<Mcolor, size_t> simple_cycle_count; 		// cycle of simple vertices
  std::map<Mcolor, size_t> special_cycle_count; 	// cycle of simple vertices and oo, of even length

  //chromosome
  std::vector<size_t> liniar_chr; 				
  std::vector<size_t> circular_chr; 				
};
#endif
