#ifndef STAGE_MANAGER_HPP
#define STAGE_MANAGER_HPP

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

	void add_prestage(algo::AbsStage<graph_pack_t>* stage) { 
		prestages.push_back(std::unique_ptr<algo::AbsStage<graph_pack_t> >(stage));
    prestages.back()->m_parent = this;
	}

  void add_stage(algo::AbsStage<graph_pack_t>* stage) {
    main_stages.push_back(std::unique_ptr<algo::AbsStage<graph_pack_t> >(stage));
    main_stages.back()->m_parent = this;
  }

  void add_poststage(algo::AbsStage<graph_pack_t>* stage) { 
		poststages.push_back(std::unique_ptr<algo::AbsStage<graph_pack_t> >(stage));
    poststages.back()->m_parent = this;
	}

	void add_prestage(std::initializer_list<AbsStage<graph_pack_t>* > stages) {
    for (auto it = stages.begin(); it = stages.end(); ++it) {
    	add_prestage(*it);
    }
  }

  void add_stage(std::initializer_list<AbsStage<graph_pack_t>* > stages) {
    for (auto it = stages.begin(); it = stages.end(); ++it) {
    	add_stage(*it);
    }
  }

	void add_poststage(std::initializer_list<AbsStage<graph_pack_t>* > stages) {
    for (auto it = stages.begin(); it = stages.end(); ++it) {
    	add_poststage(*it);
    }
  }

  void run(graph_pack_t & graph_pack);

  DebugPolicy const & get_debug_policy() const {
    return debug_policy;
  }

private:
  std::vector<std::unique_ptr<AbsStage<graph_pack_t> > > prestages;
  std::vector<std::unique_ptr<AbsStage<graph_pack_t> > > main_stages;
  std::vector<std::unique_ptr<AbsStage<graph_pack_t> > > poststages;

  DebugPolicy debug_policy;
  DECL_LOGGER("StageManager");
  size_t number_of_rounds; 
};

template<class graph_pack_t> 
void StageManager<graph_pack_t>::run(graph_pack_t& graph_pack) {
	if (debug_policy.is_save_debug) { 
		;//save stage 0 
	}

  bool isChanged = true; 
  while (isChanged) {
    isChanged = false; 

    for (size_t rounds = 1; rounds <= number_of_rounds && !isChanged; ++rounds) {  
      INFO("Start work in rounds " << rounds);
      graph_pack.update_number_of_splits(rounds);
      for (auto it_stage = main_stages.begin(); it_stage != main_stages.end() && !isChanged; ++it_stage) {        
        AbsStage<graph_pack_t> * stage = it_stage->get();
        if (rounds <= stage->get_max_round()) {        
          INFO("STAGE == " << stage->name());  
          isChanged = stage->run(graph_pack);
 
          if (isChanged && debug_policy.is_save_debug) { 
            ;//stage->save(g, saves_policy_.save_to_);
          }
        }
      }
    } 

    if (!isChanged) {
      graph_pack.update_number_of_splits(3);
      for (auto it_stage = poststages.begin(); it_stage != poststages.end() && !isChanged; ++it_stage) {        
        AbsStage<graph_pack_t> * stage = it_stage->get();
        INFO("POSTSTAGE == " << stage->name());    
        isChanged = stage->run(graph_pack);
      }
    }
  } 

  if (debug_policy.is_save_debug) { 
		;//save last graph and statistics  
	}
}

}

#endif