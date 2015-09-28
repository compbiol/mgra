#ifndef STAGE_MANAGER_HPP
#define STAGE_MANAGER_HPP

#include "writer/txt_stat.hpp" 
#include "writer/dot_graph.hpp" 

namespace algo { 

template<class graph_pack_t> 
struct StageManager {  
  struct DebugPolicy {
    bool is_save_debug;
    std::string path_to_save;

    DebugPolicy()
    : is_save_debug(false)
    , path_to_save("") 
    {
    }

    DebugPolicy(bool is_debug, std::string const & save_to)
    : is_save_debug(is_debug)
    , path_to_save(save_to) 
    {
    }
  };

  StageManager(size_t rounds = 1, DebugPolicy policy = DebugPolicy())
  : debug_policy(policy) 
  , number_of_rounds(rounds)
  {
    assert(rounds > 0 && rounds <= 3);
  }

  void add_stage(kind_stage stage);

  void add_stage(AbsStage<graph_pack_t>* stage) {
    if (stage->get_stage_type() == pre_stage_t) { 
      prestages.push_back(std::unique_ptr<algo::AbsStage<graph_pack_t> >(stage));
      prestages.back()->m_parent = this;
    } else if (stage->get_stage_type() == round_stage_t) {
      main_stages.push_back(std::unique_ptr<algo::AbsStage<graph_pack_t> >(stage));
      main_stages.back()->m_parent = this;
    } else if (stage->get_stage_type() == post_stage_t) {
      poststages.push_back(std::unique_ptr<algo::AbsStage<graph_pack_t> >(stage));
      poststages.back()->m_parent = this;
    } 
  }

  void add_stage(std::initializer_list<AbsStage<graph_pack_t>* > stages) {
    for (auto it = stages.begin(); it = stages.end(); ++it) {
    	add_stage(*it);
    }
  }

  bool run(graph_pack_t & graph_pack);

  DebugPolicy const & get_debug_policy() const {
    return debug_policy;
  }

private:
  std::vector<std::unique_ptr<AbsStage<graph_pack_t> > > prestages;
  std::vector<std::unique_ptr<AbsStage<graph_pack_t> > > main_stages;
  std::vector<std::unique_ptr<AbsStage<graph_pack_t> > > poststages;

  DebugPolicy debug_policy;
  size_t number_of_rounds; 

private: 
  DECL_LOGGER("StageManager");
};

template<class graph_pack_t> 
bool StageManager<graph_pack_t>::run(graph_pack_t& graph_pack) {
  writer::TXT_statistics<typename graph_pack_t::Statistics> debug_stat; 
  writer::GraphDot<graph_pack_t> debug_dots;

  size_t num_stage = 0;
  bool isChanged = true; 
  auto update_lambda = [&] () -> void { 
    if (isChanged && debug_policy.is_save_debug) { 
      graph_pack.update_graph_statistics();
      debug_stat.print_all_statistics(num_stage, graph_pack.stats);
      debug_dots.save_bp_graph(graph_pack, num_stage);
      ++num_stage;
    }
  };
	
  if (debug_policy.is_save_debug) { 
    debug_stat.open(debug_policy.path_to_save);
    debug_dots.open(debug_policy.path_to_save);
    debug_dots.save_subtrees();
    update_lambda();
	}

  while (isChanged) {
    isChanged = false; 

    if (!isChanged) {
      graph_pack.update_number_of_splits(1);
      for (auto it_stage = prestages.begin(); it_stage != prestages.end() && !isChanged; ++it_stage) {        
        AbsStage<graph_pack_t> * stage = it_stage->get();
        INFO("PRESTAGE == " << stage->name());    
        isChanged = stage->run(graph_pack);
        update_lambda();
      }
    }

    for (size_t rounds = 1; rounds <= number_of_rounds && !isChanged; ++rounds) {  
      INFO("Start work in rounds " << rounds);
      graph_pack.update_number_of_splits(rounds);
      
      for (auto it_stage = main_stages.begin(); it_stage != main_stages.end() && !isChanged; ++it_stage) {        
        AbsStage<graph_pack_t> * stage = it_stage->get();
        if (rounds <= stage->get_max_round()) {        
          INFO("STAGE == " << stage->name());  
          isChanged = stage->run(graph_pack); 
          update_lambda();
        }
      }
    } 

    if (!isChanged) {
      graph_pack.update_number_of_splits(3);
      for (auto it_stage = poststages.begin(); it_stage != poststages.end() && !isChanged; ++it_stage) {        
        AbsStage<graph_pack_t> * stage = it_stage->get();
        INFO("POSTSTAGE == " << stage->name());    
        isChanged = stage->run(graph_pack);
        update_lambda();
      }
    }
  } 

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

}
   
#endif