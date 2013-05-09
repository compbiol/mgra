#ifndef STAGE1_H_
#define STAGE1_H_

#include <list>
#include <string>
#include <set>

#include "mpbgraph.h"
#include "2break.h"


typedef std::list<vertex_t> path_t;

size_t process_simple_path(path_t path, MBGraph& graph);

vertex_t find_simple_path(path_t& path, const MBGraph& graph, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, bool is_next); 

// Stage 1: loop over vertices    
bool stage1(MBGraph& graph); 


#endif
