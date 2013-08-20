#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "mbgraph_history.h"
#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
	Algorithm(const std::shared_ptr<graph_t>& gr) 
	: graph(gr) 
	, canformQoo(true)
	, split_bad_colors(false)
        , max_size_component(25)
	, write_stats("stats.txt") {  
	} 

	void convert_to_identity_bgraph(const ProblemInstance<Mcolor>& cfg);

        std::map<Mcolor, std::set<arc_t> > get_removing_edges() const; 
private: 
	//Stage 1: loop over vertices  
	bool stage1(); 
	vertex_t find_simple_path(path_t& path, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, bool is_next); 	
	size_t process_simple_path(path_t& path);	

	//Stage 2: process H-subgraph
	bool stage2();
	bool canformQ(const vertex_t& x, const Mcolor& Q) const;
	bool is_mobil_edge(const vertex_t& y, const Mularcs<Mcolor>& mularcs_x, const Mularcs<Mcolor>& mularcs_y) const;

	//Stage 3: process insertion/deletion events
	bool stage3();
	//void remove_postponed_deletions();
	
	//Stage 4: process tandem duplication events
	bool stage4_td(); 
	bool stage4_rtd(); 

	//Stage 5: process components, process 4-cycles
	bool stage5_1(); 
	bool stage5_2();

	//Stage 10: convert duplication to tandem duplication   
	bool stage4_conv_to_td();

	//Stage 6: process insertion/deletion events bu splitting colors 
	//Stage 7: process H-subgraph with split bad color
	//Stage 8: process complete but non-T-consistent paths/cycles by splitting colors
	//Stage 9: process complete but non-T-consistent tandem duplication and reverse tandem duplication by splitting colors

 	//Stage 10: Force stage
	bool stage6();
	size_t calculate_cost(const vertex_t& y, const Mularcs<Mcolor>& mularcs_x, const Mularcs<Mcolor>& mulacrs_y); 
        std::set<arc_t> create_minimal_matching(const std::set<vertex_t>& vertex_set); 
	size_t process_minimal_matching(const arc_t& matching);

private: 
  std::shared_ptr<graph_t> graph; 

  bool canformQoo;  // safe choice, at later stages may change to false
  bool split_bad_colors;
  const size_t max_size_component;

  std::map<arc_t, Mcolor> postponed_deletions; 
  std::multimap<arc_t, Mcolor> insertions;

  writer::Wstats write_stats;
  writer::Wdots<graph_t, ProblemInstance<Mcolor> > write_dots; 
};

template<class graph_t>
void Algorithm<graph_t>::convert_to_identity_bgraph(const ProblemInstance<Mcolor>& cfg) {
  std::array<bool, 11> print_dots;
  print_dots.fill(true);
  bool isChanged = false;
  bool process_compl = true; 
  auto saveInfoLambda = [&](size_t stage) -> void { 
    if (print_dots[stage] && !isChanged) {
      print_dots[stage] = false; 
      Statistics<graph_t> st(graph); 
      write_stats.print_all_statistics(stage, st, cfg, *graph);
      write_dots.save_dot(*graph, cfg, stage);
    } 
  };
  saveInfoLambda(0);
   
  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    split_bad_colors = false;

#ifdef VERSION2
    if ((cfg.get_stages() >= 1) && !isChanged) { 
      std::cerr << "Stage: 1 (indel stage)" << std::endl;
      isChanged = stage3(); 
      saveInfoLambda(1);
    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      std::cerr << "Stage: 2 (indel stage)" << std::endl;
      split_bad_colors = true; 
      isChanged = stage3();
      split_bad_colors = false; 
      saveInfoLambda(2);
    }
   
    if ((cfg.get_stages() >= 3) && !isChanged) {
      std::cerr << "Stage: 3" << std::endl;
      isChanged = stage1();	
      saveInfoLambda(3);
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      std::cerr << "Stage: 4" << std::endl;
      isChanged = stage2();
      saveInfoLambda(4);
    }
   
    if ((cfg.get_stages() >= 5) && !isChanged) { // STAGE 4, somewhat unreliable
      std::cerr << "Stage: 5" << std::endl;

      isChanged = stage5_1(); // cut the graph into connected components
  
      if (!isChanged) { 
         isChanged = stage5_2(); // process 4-cycles
      }       

      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible
      }   

      saveInfoLambda(5);
    }
    
    split_bad_colors = true; 
     
    if ((cfg.get_stages() >= 6) && !isChanged) {
      std::cerr << "Stage: 6" << std::endl;
      isChanged = stage2();
      saveInfoLambda(6);
    }
   
    if ((cfg.get_stages() >= 7) && !isChanged) {
      std::cerr << "Stage: 7" << std::endl;
      isChanged = stage1();
      saveInfoLambda(7);
    }
   
    if ((cfg.get_stages() >= 8) && !isChanged) {
      std::cerr << "Stage: 8" << std::endl;
      isChanged = stage6();
      saveInfoLambda(8);
    }

    /*if ((cfg.get_stages() >= 4) && !isChanged) { 
      std::cerr << "Stage: 4 (tandem duplication stage)" << std::endl;

      isChanged = stage4_rtd(); 
      if (!isChanged) { 
	isChanged = stage4_td(); 
      }	
      if (!isChanged) { 
        isChanged = stage4_conv_to_td(); 
      }  
      saveInfoLambda(4);
    }*/


    /*if ((cfg.get_stages() >= 9) && !isChanged) {
      std::cerr << "Stage: 9" << std::endl;

      isChanged = stage4_td();
      if (!isChanged) { 
	isChanged = stage4_rtd(); 
      }	
      if (!isChanged) { 
        isChanged = stage4_conv_to_td(); 
      }
      saveInfoLambda(9);
    }*/

#else
   if ((cfg.get_stages() >= 1) && !isChanged) {
      //std::cerr << "Stage: 1" << std::endl;
      isChanged = stage1();	
      saveInfoLambda(1);
    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      //std::cerr << "Stage: 2" << std::endl;
      isChanged = stage2();
      saveInfoLambda(2);
    }

    if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 3, somewhat unreliable
      //std::cerr << "Stage: 3" << std::endl;
      isChanged = stage5_1(); // cut the graph into connected components
      
      if (!isChanged) { 
         isChanged = stage5_2(); // process 4-cycles
      }  
     
      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible
      }    

      saveInfoLambda(3);
    }

    split_bad_colors = true; 

    if ((cfg.get_stages() >= 4) && !isChanged) {
      //std::cerr << "Stage: 4" << std::endl;
      isChanged = stage2();
      saveInfoLambda(4);
    }

    if (process_compl && !cfg.get_completion().empty() && !isChanged) {     
      //std::cerr << "Manual Completion Stage" << std::endl;

      auto completion = cfg.get_completion();
      for(auto il = completion.begin(); il != completion.end(); ++il) {
	graph->apply_two_break(*il);
      }

      process_compl = false;
      isChanged = true;
    }
#endif
  }	

  write_dots.save_dot(*graph, cfg, 99);
  write_stats.print_history_statistics(*graph, get_removing_edges());
 
#ifdef VERSION3
  write_dots.save_components(*graph, cfg, 5);
#endif
}          

template<class graph_t>
std::map<Mcolor, std::set<arc_t> > Algorithm<graph_t>::get_removing_edges() const { 
  std::map<Mcolor, std::set<arc_t> > answer;
  for (const auto &edge: insertions) { 
    answer[edge.second].insert(edge.first);

    /*auto colors = graph->split_color(graph->get_complement_color(edge.second), false);
    for (const auto &col : colors) {  
      answer[col].insert(edge.first);
    }*/
  } 
  for (const auto &edge: postponed_deletions) {
    auto colors = graph->split_color(edge.second, false);
    for (const auto &col : colors) {  
      answer[col].insert(edge.first);
    }
    colors = graph->split_color(graph->get_complement_color(edge.second), false);
    for (const auto &col : colors) {  
      answer[col].insert(edge.first);
    } 
  } 
  return answer;
} 

#include "Stage1.h" 
#include "Stage2.h"
#include "Stage3.h"
#include "Stage4.h"
#include "Stage5.h"
#include "Stage6.h"

//#include "Stage99.h"

#endif
