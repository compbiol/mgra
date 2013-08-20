#ifndef WSTATS_H_
#define WSTATS_H_

#include "mbgraph_history.h"
#include "estimate.h"
#include "pconf.h"
#include "reader.h"

namespace writer { 

struct Wstats { 
  Wstats(std::string name_file);
  void print_all_statistics(int stage, Statistics<mbgraph_with_history<Mcolor> >& info, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph); 
  void print_history_statistics(const mbgraph_with_history<Mcolor>& graph, const std::map<Mcolor, std::set<arc_t> >& bad_edges);

  //void print_postponed_deletion_statistics(const std::map<arc_t, Mcolor>& postponed_deletions);
  //void print_bad_complete_edges(const mbgraph_with_history<Mcolor>& graph, const std::multimap<arc_t, Mcolor>& insertions);
private:
  void print_vertex_statistics(const std::vector<size_t>& answer);
  void print_complete_edges(const mbgraph_with_history<Mcolor>& graph); 
  void print_connected_components(const mbgraph_with_history<Mcolor>& graph);
  void print_rear_characters(const std::vector<std::string>& info);
  void print_indel_statistics(const mbgraph_with_history<Mcolor>& graph, const std::map<size_t, std::pair<Mcolor, Mcolor> >& indel); 
  void print_fair_edges(const mbgraph_with_history<Mcolor>& graph, Statistics<mbgraph_with_history<Mcolor>>& info);
	
  //void print_estimated_dist(size_t stage, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph);
private: 
  void print_start_table(size_t count_column);
  void print_close_table(bool flag = true);
private: 
  std::ofstream ofstat;
};

} 
#endif
