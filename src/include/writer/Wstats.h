#ifndef WSTATS_H_
#define WSTATS_H_

#include <string>
#include <fstream>

#include "graph_colors.h"
#include "mpbgraph.h"
#include "estimate.h"
#include "gen_dist.h"
#include "pconf.h"
#include "reader.h"

typedef std::string vertex_t;

namespace writer { 
struct Wstats { 
	Wstats(std::string name_file);
	
	inline void println(const std::string& str) { 
		ofstat << str << std::endl;
	} 

	void print_all_statistics(int stage, Statistics<MBGraph>& info, const ProblemInstance& cfg, const MBGraph& graph, const Graph_with_colors<Mcolor>& colors); 
	void print_fair_edges(const MBGraph& MBG, const Graph_with_colors<Mcolor>& colors, Statistics<MBGraph>& info);
	void histStat(const MBGraph& graph); 

	const size_t write_parametres;
private:
	void print_complete_edges(const MBGraph& MBG, const Graph_with_colors<Mcolor>& colors); 
	void print_connected_components(const MBGraph& MBG);
	void print_rear_characters(const std::vector<std::string>& info);
	void print_not_compl_characters(const std::vector<std::string>& info);
	void print_estimated_dist(size_t stage, const ProblemInstance& cfg, const MBGraph& MBG);
private: 
	void print_start_table(size_t count_column);
	void print_close_table(bool flag = true);
private: 
	std::ofstream ofstat;
};
} 
#endif
