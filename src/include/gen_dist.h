#ifndef GEN_DIST_H_
#define GEN_DIST_H_

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <vector>

#include "utility/sym_multi_hashmap.h"
//#include "utility/sym_map.h"

using namespace std; 

typedef std::string orf_t;

typedef std::string vertex_t;
typedef std::set<std::string> edge_t;
typedef std::pair<vertex_t, vertex_t> arc_t;
typedef sym_multi_hashmap<vertex_t> partgraph_t;

extern std::map<edge_t, double> BEC; // breakpoint reuse on black edges
extern std::map<edge_t, double> GEC; // breakpoint reuse on gray edges, 
extern std::map<orf_t, double> BPR; // and endpoints
extern std::set<orf_t> USE; //used verticed 
extern std::set<orf_t> USE2; //and used at least 2 times

std::vector<size_t> genome_dist(const partgraph_t& BE, const partgraph_t& GE, const partgraph_t& OE, const bool closed = false);

#endif

