#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include <list>
#include <string>

#include "2break.h"
#include "mpbgraph.h"

#include "writer/Wstats.h"
#include "writer/Wdots.h"

typedef std::list<vertex_t> path_t;

template<class graph_t>
struct Algorithm { 
	Algorithm(graph_t& gr) 
	: graph(gr) 
	, canformQoo(true)
	, split_bad_colors(false)
	, write_stats("stats.txt") {  
	} 

	void main_algorithm(const ProblemInstance& cfg);
	
	graph_t get_graph() { 	
		return graph;
	} 
private: 
	//Stage 1: loop over vertices  
	bool stage1(); 
	vertex_t find_simple_path(path_t& path, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, bool is_next); 	
	size_t process_simple_path(path_t& path);	

	//Stage 2: process H-subgraph
	bool stage2();
	bool canformQ(const std::string& x, const Mcolor& Q) const;

	//Stage 3: process components, process 4-cycles
	bool stage3_1(); 
	bool stage3_2();

	//Stage 4: process H-subgraph with split bad color

	//Not uses stage: 
	bool cut_free_ends(); 
	bool find_reliable_path(); 
private:
	//Save information
	void save_information(size_t stage, const ProblemInstance& cfg);
private: 
	graph_t graph; 

	bool canformQoo;  // safe choice, at later stages may change to false
	bool split_bad_colors;

	writer::Wstats write_stats;
	writer::Wdots write_dots; 
};

template<class graph_t>
void Algorithm<graph_t>::main_algorithm(const ProblemInstance& cfg) {
  save_information(0, cfg);
  std::array<bool, 5> print_dots;
  print_dots.fill(true);
  bool process_compl = true; 
  bool isChanged = true;

  while(isChanged) {
    isChanged = false; 

    if ((cfg.get_stages() >= 1) && !isChanged) {
      outlog << "Stage: 1" << std::endl;
      isChanged = stage1();	

      if (print_dots[1] && !isChanged) {
	print_dots[1] = false;    	
	save_information(1, cfg);		
      }

    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      outlog << "Stage: 2" << std::endl;
      isChanged = stage2();

      if (print_dots[2] && !isChanged) {
	print_dots[2] = false;
	save_information(2, cfg);		    
      }
    }

    if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 3, somewhat unreliable
      outlog << "Stage: 3" << std::endl;

      isChanged = stage3_1(); // cut the graph into connected components
      
      if (!isChanged) { 
         isChanged = stage3_2(); // process 4-cycles
      } 
      
      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible
      }    

      if (print_dots[3] && !isChanged) {
	print_dots[3] = false;
	save_information(3, cfg);
      }
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      outlog << "Stage: 4" << std::endl;

      split_bad_colors = true; 
      isChanged = stage2();
      split_bad_colors = false;

      if (print_dots[4] && !isChanged) {
	print_dots[4] = false;
	save_information(4, cfg);
      }
    }

#ifndef VERSION2
    if (process_compl && !cfg.get_completion().empty() && !isChanged) {     
      outlog << "Manual Completion Stage" << std::endl;

      auto completion = cfg.get_completion();
      for(auto il = completion.begin(); il != completion.end(); ++il) {
	TwoBreak t;
	t.OldArc[0].first = (*il)[0];
	t.OldArc[0].second = (*il)[1];
	t.OldArc[1].first = (*il)[2];
	t.OldArc[1].second = (*il)[3];
	t.MultiColor = genome_match::name_to_mcolor((*il)[4]);

	t.apply(graph, true);
      }

      Statistics<graph_t> st(graph); 
      graph.colors.update_complement_color(st.get_new_color());

      process_compl = false;
      isChanged = true;
    } 
#endif

  }	

  write_dots.save_dot(graph, cfg, 99);

#ifndef VERSION2
  Statistics<graph_t> st(graph);
  write_stats.print_fair_edges(graph, st);
#else 
  write_dots.save_components(graph, cfg, 5);
#endif

  write_stats.histStat();
}  

template<class graph_t>
void Algorithm<graph_t>::save_information(size_t stage, const ProblemInstance& cfg) { 
  Statistics<MBGraph> st(graph); 

  graph.colors.update_complement_color(st.get_new_color());

  st.count_other();   
  auto p = st.get_compl_stat();
  write_stats.print_all_statistics(stage, st, cfg, graph);
  write_dots.save_dot(graph, cfg, stage);
} 

#include "Stage1.h" 
#include "Stage2.h"
#include "Stage3.h"
#include "Stage99.h"

#endif
