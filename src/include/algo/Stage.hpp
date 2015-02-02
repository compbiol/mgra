#ifndef ABSSTAGE_HPP
#define ABSSTAGE_HPP

#include "graph/graph_pack.hpp"

namespace algo { 

template<class graph_pack_t> 
struct StageManager; 

/*
 *
 */
template<class graph_pack_t> 
struct AbsStage {  
	AbsStage(std::string const & name, std::string const & id, size_t max_round = 0)
	: m_parent(nullptr) 
	, max_avaliable_round(max_round) 
	, m_name(name)
	, m_id(id)
	{
	}

	AbsStage(AbsStage<graph_pack_t> && stage) = default;
	AbsStage(AbsStage<graph_pack_t> const & stage) = delete;
	AbsStage<graph_pack_t>& operator=(AbsStage<graph_pack_t> && stage) = default;
  AbsStage<graph_pack_t>& operator=(AbsStage<graph_pack_t> const & stage) = delete;

  std::string const & name() const { 
  	return m_name; 
  }

  std::string const & id() const { 
  	return m_id; 
  }

  size_t get_max_round() const { 
  	return max_avaliable_round;
  }

  virtual bool run(graph_pack_t & graph) = 0;

protected:
  StageManager<graph_pack_t>* m_parent;
	friend struct StageManager<graph_pack_t>;

private:
	size_t max_avaliable_round; 
  std::string m_name;
  std::string m_id;
};

/*
 *
 */
template<class graph_pack_t>
struct ChangeCanformInfinity : public AbsStage<graph_pack_t> {

	ChangeCanformInfinity(size_t max_round = 3) 
	: AbsStage<graph_pack_t>("Change canformQoo", "ChangeCanformInf", max_round)
	{
	}

	bool run(graph_pack_t & graph_pack) { 
		if (graph_pack.get_canformQoo()) { 
      graph_pack.set_canformQoo(false);
      return true; 
		} else { 
			return false;
		} 
	}

private: 
	DECL_LOGGER("ChangeCanformInf")	
};

/*
 *
 */
template<class graph_pack_t>
struct ProcessComplection : public AbsStage<graph_pack_t> {
  using transform_t = typename graph_pack_t::transform_t;

	ProcessComplection(transform_t const & complection, size_t max_round = 3) 
	: AbsStage<graph_pack_t>("Change canformQoo", "ProcessComplection", max_round)
	, is_process(false)
	, m_complection(complection)
	{
	}

	bool run(graph_pack_t & graph_pack) { 
		if (!m_complection.empty() && !is_process) { 
      is_process = true;
      for (auto const & twobreak : m_complection) {
        graph_pack.apply(twobreak);
      }
      return true;
		} else { 
			return false;
		} 
	}

private: 
	bool is_process; 
	transform_t const & m_complection;
	DECL_LOGGER("ProcessComplection")	
};

/*
 *
 */
template<class graph_pack_t> 
void split_by_mobile_property(graph_pack_t & graph_pack,
      vertex_t const & v, 
      typename graph_pack_t::mularcs_t const & mularcs, 
      std::set<typename graph_pack_t::arc_t>& mobiles, 
      std::set<typename graph_pack_t::arc_t>& non_mobiles)
{ 
  for (auto const & arc : mularcs) { 
    if (graph_pack.is_mobility_edge(v, arc.second, arc.first)) { 
      mobiles.insert(arc); 
    } else { 
      non_mobiles.insert(arc);
    }
  }   
}

} 

#include "StageManager.hpp"

#endif