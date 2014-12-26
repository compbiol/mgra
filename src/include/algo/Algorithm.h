#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "graph/breakpoint_graph.hpp"
#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
  Algorithm(std::shared_ptr<graph_t> const & gr) 
  : graph(gr) 
  , write_dots()
  {
  } 

  void init_writers(std::string const & work_dir, std::string const & graphname, bool debug) { 
    write_stats.open(work_dir, "stats.txt");
    write_dots.init(work_dir, graphname, debug);
  }

  void convert_to_identity_bgraph();

  void print_median_graphs();   
private: 
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  
  struct Stage;

public:
  struct Balance; 	

private:
  struct ProcessSimplePath;
  
  struct ProcessFairEdge; 

  struct ProcessClone;

  struct ProcessTwoBreakAndClone;

  struct IncreaseNumberComponents;

  struct BruteForce;

  struct ProcessWithBlossomV;

private: 
  std::shared_ptr<graph_t> graph; 

  writer::Wstats write_stats;
  writer::Wdots<graph_t, main_config<mcolor_t> > write_dots;
};

template<class graph_t>
void Algorithm<graph_t>::convert_to_identity_bgraph() {
  std::unordered_set<size_t> print_dots;
  size_t stage = 0; 
  bool isChanged = false;
  bool process_compl = true; 

  std::vector<std::shared_ptr<Stage> > algorithm(6);
  algorithm[0] = std::shared_ptr<Stage>(new Balance(graph));
  algorithm[1] = std::shared_ptr<Stage>(new ProcessSimplePath(graph));
  algorithm[2] = std::shared_ptr<Stage>(new ProcessFairEdge(graph));
  algorithm[3] = std::shared_ptr<Stage>(new ProcessClone(graph));
  algorithm[4] = std::shared_ptr<Stage>(new IncreaseNumberComponents(graph));

  auto const saveInfoLambda = [&](size_t st) -> void { 
    if ((print_dots.count(st) == 0) && !isChanged) {
      print_dots.insert(st);
      Statistics<graph_t> stat(graph);
      write_stats.print_all_statistics(st, stat, *graph);
      write_dots.save_dot(*graph, st);
    } 
  };

  saveInfoLambda(stage++);
  write_dots.write_legend_dot();
  
  if (cfg::get().stages >= 1) { 
    graph->update_number_of_splits(3);  
    isChanged = algorithm[0]->do_action();
    saveInfoLambda(stage++);
  }

  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    stage = 2;

    for (size_t i = 1; i <= cfg::get().rounds && !isChanged; ++i) {   
      INFO("Start work in rounds " << i);
      graph->update_number_of_splits(i);

      for (size_t j = 1; j <= cfg::get().stages && !isChanged; ++j) {
        isChanged = algorithm[j]->do_action();
        saveInfoLambda(stage++);
      } 
    } 

    if (graph->get_canformQoo() && !isChanged) { 
      INFO("Change canformQoo")
      graph->set_canformQoo(false);
      isChanged = true;
    }

    if (process_compl && !cfg::get().completion.empty() && !isChanged) {     
      INFO("Start process manual completion stage")
      for(auto il = cfg::get().completion.cbegin(); il != cfg::get().completion.cend(); ++il) {
        graph->apply(*il);
      }

      process_compl = false;
      isChanged = true;
    }  

    if ((cfg::get().is_bruteforce) && !isChanged) {
      ProcessWithBlossomV bruteforce_action(graph);
      graph->update_number_of_splits(3);
      isChanged = bruteforce_action.do_action();
      saveInfoLambda(stage++);
    }
  }	

  graph->update_number_of_splits(3);
  write_dots.save_final_dot(*graph);
}          

template<class graph_t>
void Algorithm<graph_t>::print_median_graphs() { 
  graph->update_number_of_splits(3);
  auto medians = this->graph->get_medians_colors(); 

  for (auto const & colors : medians) { 
    std::vector<std::pair<typename graph_t::edge_t, mcolor_t> > edges;
    std::unordered_set<vertex_t> marks; 
    std::cerr << edges.size() << " " << marks.size() << std::endl;
    for(auto const & v : *graph) { 
      if (marks.count(v) == 0) { 
        auto mularcs = this->graph->get_all_adjacent_multiedges_with_info(v);     
        for(auto const & arc : mularcs) { 
          //std::cerr << arc.first << " " << cfg::get().mcolor_to_name(arc.second) << std::endl;
          if (colors.count(arc.second) != 0 && marks.count(arc.first) == 0) { 
            edges.push_back(std::make_pair(typename graph_t::edge_t(v, arc.first), arc.second));
          }
        }
        marks.insert(v);
      }  
    }

    std::string name = ""; 
    for(auto const & color : colors) { 
      name += (cfg::get().mcolor_to_name(color) + "_");
    }
    name += ".dot";
    
    std::cerr << "Save " << name << std::endl;
    write_dots.save_median(name, edges); 

  }
}

#include "algo/Stage.hpp"
#include "algo/BalanceStage.hpp"
#include "algo/SimplePathStage.hpp" 
#include "algo/FairEdgeStage.hpp"
#include "algo/CloneStage.hpp"
#include "algo/BothStage.hpp"
#include "algo/IncreaseComponents.hpp"

#include "algo/BruteForce.hpp"
#include "algo/BlossomVStage.hpp"

#endif
