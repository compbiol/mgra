#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "graph/graph_pack.hpp"

#include "algo/Stage.hpp"

//#include "algo/recover_tree/RecoveredTree.hpp"

#include "algo/main/Balance.hpp"
#include "algo/main/SimplePath.hpp"
#include "algo/main/FourCycles.hpp"
#include "algo/main/FairEdge.hpp"
#include "algo/main/Clone.hpp"
#include "algo/main/IncreaseComponents.hpp"
#include "algo/main/FairEdgeAndClone.hpp"

#include "algo/main/BlossomVWrapper.hpp"
#include "algo/main/BruteForce.hpp"

//#include "algo/linearization/Linearization.hpp"

namespace algo { 

template<class graph_pack_t> 
void StageManager<graph_pack_t>::add_stage(kind_stage stage) { 
  if (stage == balance_k) { 
    add_stage(new Balance<graph_pack_t>(1));
  } else if (stage == simple_path_k) {
    add_stage(new ProcessSimplePath<graph_pack_t>(3));
  } else if (stage == four_cycles_k) {
    add_stage(new ProcessFourCycles<graph_pack_t>(1));
  } else if (stage == fair_edge_k) {
    add_stage(new ProcessFairEdge<graph_pack_t>(3));
  } else if (stage == clone_k) {
    add_stage(new ProcessClone<graph_pack_t>(3));
  } else if (stage == fair_clone_edge_k) {
    add_stage(new ProcessTwoBreakAndClone<graph_pack_t>(3));
  } else if (stage == components_k) {
    add_stage(new IncreaseNumberComponents<graph_pack_t>(1));
  } else if (stage == change_canform_k) {
    add_stage(new ChangeCanformInfinity<graph_pack_t>());
  } else if (stage == bruteforce_k) { 
    add_stage(new BruteForce<graph_pack_t>(cfg::get().size_component_in_bruteforce));
  } else if (stage == blossomv_k) {
    add_stage(new ProcessWithBlossomV<graph_pack_t>());
  } else if (stage == completion_k) {
    add_stage(new ProcessComplection<graph_pack_t>(cfg::get().completion));
  } 
}

}

/**
 * This pipeline uses in main converter, when avaliable
 * full phylogenetic tree and genomes doesn't need to scaffold
 */
template<class graph_pack_t> 
bool main_algorithm(graph_pack_t & graph_pack) {
  assert(cfg::get().how_build == default_algo); 

  INFO("Start algorithm for convert from breakpoint graph to identity breakpoint graph");

  using namespace algo; 
  StageManager<graph_pack_t> algorithm(cfg::get().rounds, {cfg::get().is_debug, cfg::get().out_path_to_debug_dir});

  for (auto const & name_stage : cfg::get().pipeline) { 
    algorithm.add_stage(name_stage); 
  }

  INFO("Run algorithms stages")
  algorithm.run(graph_pack);
    
  if (!graph_pack.graph.is_identity()) { 
    INFO("T-transformation is not complete. Cannot reconstruct genomes.")
    return false; 
  } 
  
  if (!graph_pack.is_consistency_graph()) {
    INFO("We have problem with edges, corresponding postponed deletions.")
    INFO("If you have indentity breakpoint graph after stages, please contact us.")
    return false;
  } 

  INFO("Start to replace cloning to 2-breaks")
  graph_pack.history.change_history();
  INFO("Finish to replace cloning to 2-breaks")
   
  return true;
}

template<class graph_pack_t> 
bool wgd_algorithm(graph_pack_t & graph_pack) {
  assert(cfg::get().how_build == default_algo); 

  INFO("Start wgd algorithm for convert from breakpoint graph to identity breakpoint graph");

  using namespace algo; 
  StageManager<graph_pack_t> algorithm(cfg::get().rounds, {cfg::get().is_debug, cfg::get().out_path_to_debug_dir});

  INFO("Run algorithms stages")
  algorithm.run(graph_pack);
    
  if (!graph_pack.graph.is_identity()) { 
    INFO("T-transformation is not complete. Cannot reconstruct genomes.")
    return false; 
  } 
  
  if (!graph_pack.is_consistency_graph()) {
    INFO("We have problem with edges, corresponding postponed deletions.")
    INFO("If you have indentity breakpoint graph after stages, please contact us.")
    return false;
  } 

  INFO("Start to replace cloning to 2-breaks")
  graph_pack.history.change_history();
  INFO("Finish to replace cloning to 2-breaks")
   
  return true;
}

/*
template<class graph_pack_t>
struct Algorithm { 
  using mcolor_t = typename graph_pack_t::mcolor_type; 

  Algorithm() 
  : write_dots()
  {
  } 

  void init_writers(std::string const & work_dir, std::string const & graphname, bool debug) { 
    write_stats.open(work_dir, "stats.txt");
    write_dots.init(work_dir, graphname, debug);
  }

  bool main_algorithm(graph_pack_t & graph_pack);

private: 
  writer::Wstats write_stats;
  writer::Wdots<graph_pack_t, main_config<mcolor_t> > write_dots;
};

template<class graph_pack_t>
bool Algorithm<graph_pack_t>::main_algorithm(graph_pack_t & graph_pack) {
  using namespace algo;
  size_t stage = 0; 
  bool isChanged = false;
  bool process_compl = true; 

  std::vector<std::shared_ptr<AbsStage<graph_pack_t> > > algorithm(6);
  algorithm[0] = std::shared_ptr<AbsStage<graph_pack_t> >(new Balance<graph_pack_t>(1));
  algorithm[1] = std::shared_ptr<AbsStage<graph_pack_t> >(new ProcessSimplePath<graph_pack_t>(3));
  algorithm[2] = std::shared_ptr<AbsStage<graph_pack_t> >(new ProcessFourCycles<graph_pack_t>(1));
  algorithm[3] = std::shared_ptr<AbsStage<graph_pack_t> >(new ProcessFairEdge<graph_pack_t>(3));
  algorithm[4] = std::shared_ptr<AbsStage<graph_pack_t> >(new ProcessClone<graph_pack_t>(3));
  algorithm[5] = std::shared_ptr<AbsStage<graph_pack_t> >(new IncreaseNumberComponents<graph_pack_t>(1));
  
  auto const saveInfoLambda = [&](size_t st) -> void { 
    if (isChanged) {
      //Statistics<graph_t> stat(graph);
      //write_stats.print_all_statistics(st, stat, *graph);
      write_dots.save_bp_graph(graph_pack, st);
    } 
  };

  saveInfoLambda(stage++);
  write_dots.save_subtrees();
  
  if (cfg::get().stages >= 1) { 
    graph_pack.update_number_of_splits(3);  
    isChanged = algorithm[0]->run(graph_pack);
    saveInfoLambda(stage++);
  }

  isChanged = true;
  while(isChanged) {
    isChanged = false; 
    
    for (size_t i = 1; i <= cfg::get().rounds && !isChanged; ++i) {   
      INFO("Start work in rounds " << i);
      graph_pack.update_number_of_splits(i);

      for (size_t j = 1; j <= cfg::get().stages && !isChanged; ++j) {
        isChanged = algorithm[j]->run(graph_pack);       
        //INFO("Result for working " << isChanged)
        if (isChanged) { 
          INFO("STAGE == " << algorithm[j]->name());  
        
          graph_pack.update_graph_statistics();
          saveInfoLambda(stage++);
        } 
      } 
    } 

    if (graph_pack.get_canformQoo() && !isChanged) { 
      INFO("Change canformQoo")
      graph_pack.set_canformQoo(false);
      isChanged = true;
    }

    if (process_compl && !cfg::get().completion.empty() && !isChanged) {     
      INFO("Start process manual completion stage")
      for(auto il = cfg::get().completion.cbegin(); il != cfg::get().completion.cend(); ++il) {
        graph_pack.apply(*il);
      }

      process_compl = false;
      isChanged = true;
    }  

    if ((cfg::get().is_bruteforce) && !isChanged) {
      ProcessWithBlossomV<graph_pack_t> bruteforce_action;
      graph_pack.update_number_of_splits(3);
      isChanged = bruteforce_action.run(graph_pack);
      saveInfoLambda(stage++);
    }
  }	
  
  graph_pack.update_number_of_splits(3);
  write_dots.save_final_dot(graph_pack);

  return (graph_pack.graph.is_identity());
}          
*/

#endif
