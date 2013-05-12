#ifndef STAGE1_H_
#define STAGE1_H_

#include <list>
#include <set>
#include <string>

#include "mpbgraph.h"
#include "2break.h"

typedef std::list<vertex_t> path_t;

//add graph - shared_ptr, create template_class.
struct Stage1 { 
	// Stage 1: loop over vertices   
	Stage1(MBGraph& gr)
	: graph(gr) { 
	} 

	bool stage1(/*MBGraph& graph*/); 
	
	MBGraph get_graph() {
		return graph;
 	}
private: 
	size_t process_simple_path(/*path_t& path, MBGraph& graph*/);	
	vertex_t find_simple_path(/*path_t& path, MBGraph& graph,*/ const vertex_t& prev, const vertex_t& cur, bool is_next); 
private: 
	path_t path;
	std::unordered_set<vertex_t> processed;
	MBGraph graph;
};


#endif
