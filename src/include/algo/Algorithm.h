#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include <list>
#include <string>

#include "mbgraph_history.h"

#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
	Algorithm(const std::shared_ptr<graph_t>& gr) 
	: graph(gr) 
	, canformQoo(true)
	, split_bad_colors(false)
	, write_stats("stats.txt") {  
	} 

	void convert_to_identity_bgraph(const ProblemInstance<Mcolor>& cfg);

private: 
	//Stage 1: loop over vertices  
	bool stage1(); 
	vertex_t find_simple_path(path_t& path, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, bool is_next); 	
	size_t process_simple_path(path_t& path);	

	//Stage 2: process H-subgraph
	bool stage2();
	bool canformQ(const std::string& x, const Mcolor& Q) const;
	bool is_mobil_edge(const vertex_t& y, const Mularcs<Mcolor>& Cx, const Mularcs<Mcolor>& Cy) const;

	//Stage 3: process insertion/deletion events
	bool stage3_1();
	void remove_past_bad_colors();

	//Stage 4: process tandem duplication events
	bool stage4_td(); 
	bool stage4_rtd(); 

	//Stage 5: process components, process 4-cycles
	bool stage5_1(); 
	bool stage5_2();

	//Stage 10: convert duplication to tandem duplication  //FIXME: GO TO STAGE 6 
	bool stage4_conv_to_td();

	//Stage 6: process insertion/deletion events bu splitting colors 
	//Stage 7: process H-subgraph with split bad color
	//Stage 8: process complete but non-T-consistent paths/cycles by splitting colors
	//Stage 9: process complete but non-T-consistent tandem duplication and reverse tandem duplication by splitting colors

	//Not uses stage: 
	//bool cut_free_ends(); 
	//bool find_reliable_path(); 
private:
	//Save information
	void save_information(size_t stage, const ProblemInstance<Mcolor>& cfg);
private: 
	std::shared_ptr<graph_t> graph; 

	bool canformQoo;  // safe choice, at later stages may change to false
	bool split_bad_colors;

	std::list<InsDel<Mcolor> > viewed_edges; // , Mcolor viewed edges for stage3 
	
	writer::Wstats write_stats;
	writer::Wdots<graph_t, ProblemInstance<Mcolor> > write_dots; 
};

template<class graph_t>
void Algorithm<graph_t>::convert_to_identity_bgraph(const ProblemInstance<Mcolor>& cfg) {
  save_information(0, cfg);
  std::array<bool, 11> print_dots;
  print_dots.fill(true);
  bool process_compl = true; 
  bool isChanged = true;
 
  write_dots.write_legend_dot(cfg);
  while(isChanged) {
    isChanged = false; 

#ifdef VERSION2
    if ((cfg.get_stages() >= 1) && !isChanged) {
      std::cerr << "Stage: 1" << std::endl;

      isChanged = stage1();	

      if (print_dots[1] && !isChanged) {
	print_dots[1] = false;    	
	save_information(1, cfg);		
      }

    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      std::cerr << "Stage: 2" << std::endl;

      isChanged = stage2();

      if (print_dots[2] && !isChanged) {
	print_dots[2] = false;
	save_information(2, cfg);		    
      }
    }

    if ((cfg.get_stages() >= 3) && !isChanged) { 
      std::cerr << "Stage: 3 (indel stage)" << std::endl;

      if (!isChanged) { 
	isChanged = stage3_1(); 
      }	

      if (print_dots[3] && !isChanged) {
	print_dots[3] = false;
	save_information(3, cfg);		    
      }
    }


   if ((cfg.get_stages() >= 4) && !isChanged) { 
      std::cerr << "Stage: 4 (tandem duplication stage)" << std::endl;

      isChanged = stage4_rtd(); 
      
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

      split_bad_colors = true; 
      isChanged = stage3_1();
      split_bad_colors = false;

      if (print_dots[6] && !isChanged) {
	print_dots[6] = false;
	save_information(6, cfg);
      }
    }

    /*if ((cfg.get_stages() >= 7) && !isChanged) {
      std::cerr << "Stage: 7" << std::endl;

      split_bad_colors = true; 
      isChanged = stage1();
      split_bad_colors = false;

      if (print_dots[7] && !isChanged) {
	print_dots[7] = false;
	save_information(7, cfg);
      }
    }*/
   
    if ((cfg.get_stages() >= 8) && !isChanged) {
      std::cerr << "Stage: 8" << std::endl;

      split_bad_colors = true; 
      isChanged = stage2();
      split_bad_colors = false;

      if (print_dots[8] && !isChanged) {
	print_dots[8] = false;
	save_information(8, cfg);
      }
    }
    
    /*if ((cfg.get_stages() >= 9) && !isChanged) {
      std::cerr << "Stage: 9" << std::endl;

      split_bad_colors = true; 
      isChanged = stage4_td();
      split_bad_colors = false;

      if (!isChanged) { 
	split_bad_colors = true; 
	isChanged = stage4_rtd(); 
	split_bad_colors = false;
      }	

      if (print_dots[9] && !isChanged) {
	print_dots[9] = false;
	save_information(9, cfg);
      }
    }*/


    /*if ((cfg.get_stages() >= 10) && !isChanged) { 
      std::cerr << "Stage: 10 (convert from duplication to tandem duplication)" << std::endl;

      isChanged = stage4_conv_to_td(); 
      
      if (print_dots[10] && !isChanged) {
	print_dots[10] = false;
	save_information(10, cfg);		    
      }
    }

    if ((cfg.get_stages() >= 11) && !isChanged) { 
      std::cerr << "Stage: 11 (convert from duplication to tandem duplication) less reliable" << std::endl;

      split_bad_colors = true; 
      isChanged = stage4_conv_to_td(); 
      split_bad_colors = false; 
	
      if (print_dots[11] && !isChanged) {
	print_dots[11] = false;
	save_information(11, cfg);		    
      }
    }*/
#else
   if ((cfg.get_stages() >= 1) && !isChanged) {
      //std::cerr << "Stage: 1" << std::endl;
      isChanged = stage1();	

      if (print_dots[1] && !isChanged) {
	print_dots[1] = false;    	
	save_information(1, cfg);		
      }

    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      //std::cerr << "Stage: 2" << std::endl;

      isChanged = stage2();

      if (print_dots[2] && !isChanged) {
	print_dots[2] = false;
	save_information(2, cfg);		    
      }
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

      if (print_dots[3] && !isChanged) {
	print_dots[3] = false;
	save_information(3, cfg);
      }
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      //std::cerr << "Stage: 4" << std::endl;

      split_bad_colors = true; 
      isChanged = stage2();
      split_bad_colors = false;

      if (print_dots[4] && !isChanged) {
	print_dots[4] = false;
	save_information(4, cfg);
      }
    }

#ifndef VERSION21
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
#endif

  }	

#ifdef VERSION2 
  if (!viewed_edges.empty()) {
 	remove_past_bad_colors();
  }

  if (!viewed_edges.empty()) {
	std::cerr << "WARNING, WARNING: Viewed edges, when we insert TC color, not removed. We have " << viewed_edges.size() << std::endl;
  }
#endif 

  write_dots.save_dot(*graph, cfg, 99);

#ifndef VERSION2
  Statistics<graph_t> st(graph);
  write_stats.print_fair_edges(*graph, st);
#else 
  write_dots.save_components(*graph, cfg, 5);
#endif

  write_stats.histStat(*graph);
}  

template<class graph_t>
void Algorithm<graph_t>::save_information(size_t stage, const ProblemInstance<Mcolor>& cfg) { 
  Statistics<graph_t> st(graph); 

  st.count_other();   
  auto p = st.get_compl_stat();
  write_stats.print_all_statistics(stage, st, cfg, *graph);
  write_dots.save_dot(*graph, cfg, stage);
} 

#include "Stage1.h" 
#include "Stage2.h"
#include "Stage3.h"
#include "Stage4.h"
#include "Stage5.h"

//#include "Stage99.h"

#endif
