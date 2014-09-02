#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "graph/breakpoint_graph.hpp"
#include "writer/Wstats.h"
#include "writer/Wdots.h"

template<class graph_t>
struct Algorithm { 
  Algorithm(std::shared_ptr<graph_t> const & gr, ProblemInstance<typename graph_t::mcolor_type> const & cfg) 
  : graph(gr) 
  , rounds(cfg.get_max_number_of_split())
  , max_size_component(cfg.get_size_component_in_brutforce())
  , stages(cfg.get_stages())
  , m_completion(cfg.get_completion())
  , write_dots(cfg)
  {
  } 

  void init_writers(fs::path const & work_dir, std::string const & graphname, bool debug) { 
    write_stats.open(work_dir, "stats.txt");
    write_dots.init(work_dir, graphname, debug);
  }

  void convert_to_identity_bgraph();
  
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

  size_t const rounds;
  size_t const max_size_component; 
  size_t const stages;
  std::list<twobreak_t> const m_completion;

  writer::Wstats write_stats;
  writer::Wdots<graph_t, ProblemInstance<mcolor_t> > write_dots;
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
  //algorithm[4] = std::shared_ptr<Stage>(new ProcessTwoBreakAndClone(graph));
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
  
  if (stages >= 1) { 
    graph->update_number_of_splits(3);  
    std::cerr << "Stage " << stage << ": " << algorithm[0]->get_name() << std::endl; 
    isChanged = algorithm[0]->do_action();
    saveInfoLambda(stage++);
  }

  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    stage = 2;

    for (size_t i = 1; i <= rounds && !isChanged; ++i) {   
      std::cerr << "Rounds " << i << std::endl;

      graph->update_number_of_splits(i);

      for (size_t j = 1; j < 5 && !isChanged; ++j) {
        std::cerr << "Stage " << stage << ": " << algorithm[j]->get_name() << std::endl;
        isChanged = algorithm[j]->do_action();
        saveInfoLambda(stage++);
      } 
    } 

    if (graph->get_canformQoo() && !isChanged) { 
      std::cerr << "Change canformQ2oo "  << std::endl;
      graph->set_canformQoo(false);
      isChanged = true;
    }

    if (process_compl && !m_completion.empty() && !isChanged) {     
      //std::cerr << "Manual Completion Stage" << std::endl;
      for(auto il = m_completion.cbegin(); il != m_completion.cend(); ++il) {
        graph->apply(*il);
      }

      process_compl = false;
      isChanged = true;
    }  

    if ((max_size_component != 0) && !isChanged) {
      ProcessWithBlossomV bruteforce_action(graph);
      std::cerr << "Stage " << stage << " " << bruteforce_action.get_name() << std::endl;
      graph->update_number_of_splits(3);
      isChanged = bruteforce_action.do_action();
      saveInfoLambda(stage++);
    }

  }	
       
  graph->update_number_of_splits(3);
  write_dots.save_final_dot(*graph);
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
