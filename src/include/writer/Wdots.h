#ifndef WDOTS_H_
#define WDOTS_H_

#include <string>
#include <vector>
#include <fstream>

#include "mpbgraph.h"
#include "mcolor.h"
#include "genome.h"
#include "pconf.h"


namespace writer {
struct Wdots { 
	// Save .dot file and output statistics of synteny blocks representing breakpoints
	void save_dot(const MBGraph& graph, const Graph_with_colors<Mcolor>& colors, const ProblemInstance& cfg, size_t stage);
	void save_components(const MBGraph& graph, const Graph_with_colors<Mcolor>& colors, const ProblemInstance& cfg, size_t stage);
	//void write_legend_dot(size_t size, const std::vector<std::string>& output, const ProblemInstance& cfg);
}; 
} 

#endif

