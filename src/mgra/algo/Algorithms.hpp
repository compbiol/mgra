#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "boost/optional.hpp"

#include "graph/graph_pack.hpp"
#include "algo/Stage.hpp"

#include "algo/stages/Balance.hpp"
#include "algo/stages/SimplePath.hpp"
#include "algo/stages/FourCycles.hpp"
#include "algo/stages/FairEdge.hpp"
#include "algo/stages/Clone.hpp"
#include "algo/stages/IncreaseComponents.hpp"
#include "algo/stages/FairEdgeAndClone.hpp"

#include "algo/stages/BlossomVWrapper.hpp"
#include "algo/stages/BruteForce.hpp"

#include "../linearization/algo/recovered_information.hpp"

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


/**
 * This pipeline uses in main converter, when avaliable
 * full phylogenetic tree and genomes doesn't need to scaffold
 */
template<class graph_pack_t> 
boost::optional<typename algo::RecoveredInformation<graph_pack_t>::AncestorInformation> main_algorithm(graph_pack_t & graph_pack) {
  assert(cfg::get().how_build == default_algo); 

  INFO("Start algorithm for convert from breakpoint graph to identity breakpoint graph");

  INFO("Init pipeline") 
  using namespace algo; 
  StageManager<graph_pack_t> algorithm(cfg::get().rounds, cfg::get().out_path_to_saves_dir, {cfg::get().is_debug, cfg::get().out_path_to_debug_dir});

  for (auto const &name_stage : cfg::get().pipeline) {
    algorithm.add_stage(name_stage);
  }

  INFO("Run pipeline") 
  bool is_can_reconstruct = algorithm.run(graph_pack);
  
  if (is_can_reconstruct) { 
    INFO("Get results from graphs.")
    RecoveredInformation<graph_pack_t> recover_info(graph_pack);
    if (!cfg::get().target.empty()) {
      recover_info.init_target_results();
    } else {
      recover_info.init_raw_results();
    }
    INFO("Finish get results from graphs.")
    
    return recover_info.get_results();
  } else { 
    return boost::none;
  }
}

template<class graph_pack_t> 
boost::optional<typename algo::RecoveredInformation<graph_pack_t>::AncestorInformation> wgd_algorithm(graph_pack_t & graph_pack) {
  assert(cfg::get().how_build == default_algo); 

  INFO("Start wgd algorithm for convert from breakpoint graph to identity breakpoint graph");

  using namespace algo; 
  StageManager<graph_pack_t> algorithm(cfg::get().rounds, cfg::get().out_path_to_saves_dir,  {cfg::get().is_debug, cfg::get().out_path_to_debug_dir});

  INFO("Run algorithms stages")
  algorithm.run(graph_pack);
     
  return boost::none;
}

}

#endif
