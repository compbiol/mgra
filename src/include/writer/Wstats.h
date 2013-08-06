#ifndef WSTATS_H_
#define WSTATS_H_

#include <string>
#include <fstream>

#include "mbgraph_history.h"
#include "estimate.h"
//#include "gen_dist.h"
#include "pconf.h"
#include "reader.h"

typedef std::string vertex_t;

namespace writer { 
struct Wstats { 
	Wstats(std::string name_file);
	
	inline void println(const std::string& str) { 
		ofstat << str << std::endl;
	} 

	void print_all_statistics(int stage, Statistics<mbgraph_with_history<Mcolor> >& info, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph); 
	void print_fair_edges(const mbgraph_with_history<Mcolor>& MBG, Statistics<mbgraph_with_history<Mcolor>>& info);
	void histStat(const mbgraph_with_history<Mcolor>& graph); 

	void print_duplication_statistics(const std::vector<size_t>& answer);

	const size_t write_parametres;
private:
	void print_complete_edges(const mbgraph_with_history<Mcolor>& MBG); 
	void print_connected_components(const mbgraph_with_history<Mcolor>& MBG);
	void print_rear_characters(const std::vector<std::string>& info);
	void print_not_compl_characters(const std::vector<std::string>& info);
	void print_estimated_dist(size_t stage, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& MBG);
private: 
	void print_start_table(size_t count_column);
	void print_close_table(bool flag = true);
private: 
	std::ofstream ofstat;
};
} 
#endif
