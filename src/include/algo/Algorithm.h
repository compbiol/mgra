#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include <list>
#include <string>

#include "mbgraph_history.h"

#include "writer/Wstats.h"
#include "writer/Wdots.h"

extern std::ofstream outlog;

typedef std::list<vertex_t> path_t;

template<class graph_t>
struct Algorithm { 
	Algorithm(graph_t& gr) //FIXME: get const
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

	//Stage 3: process insertion/deletion events
	bool stage3_1();
	bool stage3_2();

	//Stage 4: process tandem duplication events
	bool stage4_td(); 
	bool stage4_rtd(); 

	//Stage 5: process components, process 4-cycles
	bool stage5_1(); 
	bool stage5_2();

	//Stage 6: process insertion/deletion events less reliable
	bool stage6();

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

	std::list<std::pair<InsDel<Mcolor>, std::set<Mcolor> > > viewed_edges; // , Mcolor viewed edges for stage3_1 
	
	writer::Wstats write_stats;
	writer::Wdots write_dots; 
};

template<class graph_t>
void Algorithm<graph_t>::main_algorithm(const ProblemInstance& cfg) {
  save_information(0, cfg);
  std::array<bool, 7> print_dots;
  print_dots.fill(true);
  bool process_compl = true; 
  bool isChanged = true;

  while(isChanged) {
    isChanged = false; 

    if ((cfg.get_stages() >= 1) && !isChanged) {
#ifdef VERSION2
      std::cerr << "Stage: 1" << std::endl;
#endif
      isChanged = stage1();	

      if (print_dots[1] && !isChanged) {
	print_dots[1] = false;    	
	save_information(1, cfg);		
      }

    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
#ifdef VERSION2
      std::cerr << "Stage: 2" << std::endl;
#endif
      isChanged = stage2();

      if (print_dots[2] && !isChanged) {
	print_dots[2] = false;
	save_information(2, cfg);		    
      }
    }

#ifdef VERSION2
    if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 3
      std::cerr << "Stage: 3 (indel stage)" << std::endl;

      if (!viewed_edges.empty() && !isChanged) {
	 isChanged = stage3_2();
      }

      if (!isChanged) { 
	isChanged = stage3_1(); 
      }	

      if (print_dots[3] && !isChanged) {
	print_dots[3] = false;
	save_information(3, cfg);		    
      }
    }

   if ((cfg.get_stages() >= 4) && !isChanged) { // STAGE 4
      std::cerr << "Stage: 4 (tandem duplication stage)" << std::endl;

      if (!isChanged) { 
	isChanged = stage4_rtd(); 
      }	

      if (!isChanged) { 
	isChanged = stage4_td(); 
      }	

      if (print_dots[4] && !isChanged) {
	print_dots[4] = false;
	save_information(4, cfg);		    
      }
    }

    if ((cfg.get_stages() >= 5) && !isChanged) { // STAGE 5, somewhat unreliable
      std::cerr << "Stage: 5" << std::endl;

      isChanged = stage5_1(); // cut the graph into connected components
  
      if (!isChanged) { 
         isChanged = stage5_2(); // process 4-cycles
      }       

      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible
      }   

      if (print_dots[5] && !isChanged) {
	print_dots[5] = false;
	save_information(5, cfg);
      }
    }
    
    if ((cfg.get_stages() >= 6) && !isChanged) {
      std::cerr << "Stage: 6" << std::endl;

      isChanged = stage6();

      if (print_dots[6] && !isChanged) {
	print_dots[6] = false;
	save_information(6, cfg);
      }
    }
#else

    if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 3, somewhat unreliable
      outlog << "Stage: 3" << std::endl;

      isChanged = stage5_1(); // cut the graph into connected components
      
      if (!isChanged) { 
         isChanged = stage5_2(); // process 4-cycles
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


    if (process_compl && !cfg.get_completion().empty() && !isChanged) {     
      outlog << "Manual Completion Stage" << std::endl;

      auto completion = cfg.get_completion();
      for(auto il = completion.begin(); il != completion.end(); ++il) {
	TwoBreak<Mcolor> t((*il)[0], (*il)[1], (*il)[2], (*il)[3], genome_match::name_to_mcolor((*il)[4]));
	graph.apply_two_break(t, true);
      }

      Statistics<graph_t> st(graph); 
      graph.update_complement_color(st.get_new_color());

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
  if (!viewed_edges.empty()) {
	std::cerr << "WARNING, WARNING: Viewed edges, when we insert TC color, not removed. We have " << viewed_edges.size() << std::endl;
  }
  //write_dots.save_components(graph, cfg, 5);
#endif

  write_stats.histStat(graph);
}  

template<class graph_t>
void Algorithm<graph_t>::save_information(size_t stage, const ProblemInstance& cfg) { 
  Statistics<graph_t> st(graph); 

  graph.update_complement_color(st.get_new_color());

  st.count_other();   
  auto p = st.get_compl_stat();
  write_stats.print_all_statistics(stage, st, cfg, graph);
  write_dots.save_dot(graph, cfg, stage);
} 

#include "Stage1.h" 
#include "Stage2.h"
#include "Stage3.h"
#include "Stage4.h"
#include "Stage5.h"

#include "Stage99.h"

#endif
