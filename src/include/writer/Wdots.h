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
	//Wdots(std::string name_file);
	
	// Save .dot file and output statistics of synteny blocks representing breakpoints
	void save_dot(const MBGraph& graph, const ColorsGraph<Mcolor>& colors, const ProblemInstance& cfg, size_t stage);
	void save_components(const MBGraph& graph, const ColorsGraph<Mcolor>& colors, const ProblemInstance& cfg, size_t stage);
	void write_legend_dot(size_t size_genomes, const std::vector<std::string>& info);
	//write_legend_dot(const std::vector<Genome>& genomes, const std::vector<std::string>& output); 
}; 
} 

#endif

