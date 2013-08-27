#ifndef WSTATS_H_
#define WSTATS_H_

#include "mbgraph_history.h"
#include "estimate.h"
#include "pconf.h"
#include "reader.h"

namespace writer { 

struct Wstats { 
  typedef structure::Mcolor mcolor_t;
  typedef structure::Mularcs<mcolor_t> mularcs_t;

  Wstats(std::string name_file) 
  : ofstat(name_file)
  { 
  } 

  void print_all_statistics(size_t stage, Statistics<mbgraph_with_history<mcolor_t> >& info, const ProblemInstance<mcolor_t>& cfg, const mbgraph_with_history<mcolor_t>& graph); 
  void print_history_statistics(const mbgraph_with_history<mcolor_t>& graph, const edges_t& bad_edges);

  //void print_postponed_deletion_statistics(const std::map<arc_t, Mcolor>& postponed_deletions);
  //void print_bad_complete_edges(const mbgraph_with_history<Mcolor>& graph, const std::multimap<arc_t, Mcolor>& insertions);
private:
  inline void print_vertex_statistics(const std::array<size_t, 4>& answer) {
    ofstat << "... Duplication vertex: " << answer[0] << std::endl;
    ofstat << "... Insertion/deletion vertex: " << answer[1] << std::endl;
    ofstat << "... Count self loop: " << answer[2] << std::endl;
    ofstat << "... Colors is not one-to-one match: " << answer[3] << std::endl; 
  } 

  void print_complete_edges(const mbgraph_with_history<mcolor_t>& graph); 
  void print_connected_components(const mbgraph_with_history<mcolor_t>& graph);
  void print_rear_characters(const std::vector<std::string>& info);
  void print_indel_statistics(const std::multimap<size_t, std::tuple<mcolor_t, mcolor_t, size_t, size_t> >& indel); 

  //void print_fair_edges(const mbgraph_with_history<Mcolor>& graph, Statistics<mbgraph_with_history<Mcolor>>& info);	
  //void print_estimated_dist(size_t stage, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph);
private: 
  void print_start_table(size_t count_column) {
    ofstat << "\\begin{table}[h]" << std::endl;
    ofstat << "\\centering \\begin{tabular}{|c|";
    for(size_t i = 0; i < count_column - 1; ++i) { 
      ofstat << "c|";
    } 
    ofstat << "}" << std::endl;
    ofstat << "\\hline" << std::endl;	
  } 

  void print_close_table() { 
    ofstat << "\\hline" << std::endl;
    ofstat << "\\end{tabular}" << std::endl;
    ofstat << "\\end{table}" << std::endl;
    ofstat << std::endl;
  } 
private: 
  std::ofstream ofstat;
};

} 
#endif
