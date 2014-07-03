#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "genome_match.h" //FIXME REMOVE LATER

#include "graph/mbgraph_history.h"
#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
  Algorithm(std::shared_ptr<graph_t> const & gr, fs::path const & work_dir, size_t size_component, size_t rds, 
    std::string const & colorscheme, std::string const & graphname) 
  : graph(gr) 
  , canformQoo(true)
  , rounds(rds)
  , max_size_component(size_component)
  , write_stats(work_dir, "stats.txt") 
  , write_dots(work_dir, colorscheme, graphname)
  {
  } 

  void convert_to_identity_bgraph(ProblemInstance<typename graph_t::mcolor_t> const & cfg);
  edges_t get_bad_edges() const; 

private: 
  typedef typename graph_t::mcolor_t mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  typedef typename graph_t::fake_twobreak_t fake_twobreak_t;        
  typedef typename graph_t::insertion_t insertion_t;
  typedef typename graph_t::tandem_duplication_t tandem_duplication_t;

  //Stage 1: process insertion/deletion events. Balanced graph
  bool stage3();
	
  //Stage 2: Simple paths  
  bool stage1(); 
  size_t process_simple_path(path_t& path);	

  //Stage 3: process non-mobile edges
  bool stage2();

  bool stage22();
  bool is_mobility_edge(vertex_t const & x, vertex_t const & y) const;
  bool is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;
  size_t is_mobility_edge_score(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const;

  //Stage -: process tandem duplication events
  bool stage4_td(); 
  bool stage4_rtd(); 

  //Stage 4: process components, 
  utility::equivalence<vertex_t> split_on_components(std::map<vertex_t, std::set<arc_t> >& EC, std::map<vertex_t, std::set<arc_t> >& EI, mcolor_t const & Q);
  bool stage5_1(); 
  bool stage5_2();

  //Stage 5: Clone approach
  bool stage7();
  bool stage71();
        
  //Stage 10: convert duplication to tandem duplication   
  bool stage4_conv_to_td();

  //Stage 6: brutoforce stage
  bool stage6();
  size_t calculate_cost(const vertex_t& y, const mularcs_t& mularcs_x, const mularcs_t& mulacrs_y); 
  std::multimap<size_t, arc_t> create_minimal_matching(const std::set<vertex_t>& vertex_set); 
  //size_t process_minimal_matching(const std::set<arc_t>& matchings);
  size_t take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y);

  size_t check_postponed_deletions() const;

  //can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
  //if return false, then Q cannot be formed
  //if true - who knows... 
  bool canformQ(vertex_t const & current, mcolor_t  const & Q) const;

private: 
  std::shared_ptr<graph_t> graph; 

  bool canformQoo;  // safe choice, at later stages may change to false
  size_t const rounds;
  size_t const max_size_component; 

  edges_t insertions;
  edges_t postponed_deletions; 
  std::map<std::pair<vertex_t, mcolor_t>, vertex_t> mother_verteces;
 
  std::set<edge_t> clone_edges; 
  std::set<std::set<mcolor_t> > disjoint_subset_colors;  //FIXME

  writer::Wstats write_stats;
  writer::Wdots<graph_t, ProblemInstance<mcolor_t> > write_dots;
};

template<class graph_t>
void Algorithm<graph_t>::convert_to_identity_bgraph(ProblemInstance<mcolor_t> const & cfg) {
  std::unordered_set<size_t> print_dots;
  size_t stage = 0; 
  bool isChanged = false;
  bool process_compl = true; 

  auto const saveInfoLambda = [&](size_t st) -> void { 
    if ((print_dots.count(st) == 0) && !isChanged) {
      print_dots.insert(st);
      Statistics<graph_t> stat(graph); 
      write_stats.print_all_statistics(st, stat, *graph);
      write_dots.save_dot(*graph, cfg, st);
    } 
  };

  saveInfoLambda(stage++);
  write_dots.write_legend_dot(cfg);

#ifndef VERSION1
  if (cfg.get_stages() >= 1) { 
    std::cerr << "Stage: 1 (indel stage)" << std::endl; 
    graph->update_number_of_splits(rounds);  
    stage3();
    saveInfoLambda(stage++);
  }
#endif

  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    stage = 2; 

#ifndef VERSION1
    
    for (size_t i = 1; i <= rounds && !isChanged; ++i) {   
      std::cerr << "Rounds " << i << std::endl;

      graph->update_number_of_splits(i);
    
      if ((cfg.get_stages() >= 2) && !isChanged) {
        std::cerr << "Stage: 2 Good path " << stage << std::endl;

        isChanged = stage1();	
 
        if (!isChanged) {
          std::cerr << "Stage: 2 Non-mobile edges " << stage << std::endl;
          isChanged = stage22();
        }

        if (!isChanged) {
          std::cerr << "Stage: 71 Non-mobile edges " << stage << std::endl;
          isChanged = stage71();
        }

        saveInfoLambda(stage++);
      }

      if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 4, somewhat unreliable
        std::cerr << "Stage: 3 Split on components " << stage << std::endl;
        isChanged = stage5_1(); // cut the graph into connected components
        saveInfoLambda(stage++); 
      }

      if ((cfg.get_stages() >= 4) && !isChanged) {
        std::cerr << "Stage: 4 Clone approach " << stage << std::endl;
        isChanged = stage7();
        saveInfoLambda(stage++);
      }
    } 

    if (canformQoo && !isChanged) { 
      std::cerr << "Change canformQoo "  << std::endl;
      canformQoo = false; // more flexible
      isChanged = true;
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

    if ((max_size_component != 0) && !isChanged) {
      std::cerr << "Brute force stage: " << std::endl;
      graph->update_number_of_splits(3);
      isChanged = stage6();
      saveInfoLambda(stage++);
    }

#else

   stage = 2; 

   graph->update_number_of_splits(1);
    
    if ((cfg.get_stages() >= 1) && !isChanged) {
      std::cerr << "Stage: 1" << std::endl;
      isChanged = stage1();	
      saveInfoLambda(stage++);
    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      std::cerr << "Stage: 2" << std::endl;
      isChanged = stage2();
      saveInfoLambda(stage++);
    }

    if ((cfg.get_stages() >= 3) && !isChanged) { // STAGE 3, somewhat unreliable
      std::cerr << "Stage: 3" << std::endl;
      isChanged = stage5_1(); // cut the graph into connected components
      
      if (!isChanged) { 
         isChanged = stage5_2(); // process 4-cycles
      }  
     
      if (canformQoo && !isChanged) {
        isChanged = true;
        canformQoo = false; // more flexible
      }    

      saveInfoLambda(stage++);
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      std::cerr << "Stage: 4" << std::endl;
      graph->update_number_of_splits(3);
      isChanged = stage2();
      saveInfoLambda(stage++);
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

  graph->change_history();

  size_t bad_postponed_deletions = check_postponed_deletions();
  if (bad_postponed_deletions != 0) {
    std::cerr << "We have problem with " << bad_postponed_deletions << " edges, corresponding postponed deletions." << std::endl;
    std::cerr << "If you have indentity breakpoint graph after stages, please contact us." << std::endl;
    exit(1);
  } 

  write_stats.print_history_statistics(*graph, get_bad_edges());
}          

template<class graph_t>
edges_t Algorithm<graph_t>::get_bad_edges() const { 
  edges_t answer; 

  for (auto const & edge: insertions) { 
    answer.insert(edge.first, edge.second);
  } 

  for (auto const & edge: postponed_deletions) {
    answer.insert(edge.first, edge.second);
  }

  return answer;
} 

#include "Stage1.h" 
#include "Stage2.h"
#include "Stage3.h"
//#include "Stage4.h"
#include "Stage5.h"
#include "Stage6.h"
#include "Stage7.h"

#endif
