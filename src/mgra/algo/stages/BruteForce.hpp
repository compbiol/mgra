#ifndef BRUTE_FORCE_HPP 
#define BRUTE_FORCE_HPP

namespace algo { 

template<class graph_pack_t>
struct BruteForce : public algo::AbsStage<graph_pack_t> {
  using mcolor_t = typename graph_pack_t::mcolor_t;
  
  using edge_t = typename graph_pack_t::edge_t;  
  using arc_t = typename graph_pack_t::arc_t; 
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using twobreak_t = typename graph_pack_t::twobreak_t;

  BruteForce(size_t size_component)
  : AbsStage<graph_pack_t>("Resolve component with score function", "brute_force", post_stage_t, 3) 
  , max_size_component(size_component)
  {
  }
    
  bool run(graph_pack_t & graph_pack) override;
  
private:
  std::multimap<size_t, edge_t> create_minimal_matching(graph_pack_t const & graph_pack, 
                                                        std::set<vertex_t> const & vertex_set); 
  
private: 
  size_t max_size_component;

private:
  DECL_LOGGER("BruteForce");
};

template<class graph_pack_t>
bool BruteForce<graph_pack_t>::run(graph_pack_t & graph_pack) { 
  size_t number_rear = 0; // number of rearrangements 

  utility::equivalence<vertex_t> components = graph_pack.split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();

  for (auto const & vertex_set : classes) { 
     if ((vertex_set.second.size() <= max_size_component)) { 
      std::multimap<size_t, edge_t> matching;
      matching = create_minimal_matching(graph_pack, vertex_set.second);
      auto edge = matching.cbegin()->second;
      auto result = graph_pack.take_edge_on_color(edge.first, graph_pack.multicolors.get_complete_color(), edge.second);

      for (twobreak_t const & twobreak : result) { 
        graph_pack.apply(twobreak);
        ++number_rear;
      }      
    } 
  }

  return (number_rear != 0);
} 

template<class graph_pack_t>
std::multimap<size_t, typename graph_pack_t::edge_t> 
      BruteForce<graph_pack_t>::create_minimal_matching(graph_pack_t const & graph_pack, 
                                                  std::set<vertex_t> const & vertex_set) {
  std::map<edge_t, std::pair<size_t, mcolor_t> > weight_edges; 

  for (vertex_t const & v : vertex_set) { 
    mularcs_t mularcs = graph_pack.get_all_adjacent_multiedges(v);
    for (arc_t const & arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
      	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(graph_pack.calculate_cost(v, arc.first), arc.second)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<edge_t> current; 
  for (auto const & edge : weight_edges) { 
    bool is_have_tc = false; 
    auto edge_colors = graph_pack.split_color(edge.second.second);
    
    for (auto const & color : edge_colors) { 
      if (graph_pack.multicolors.is_T_consistent_color(color) && !graph_pack.multicolors.is_vec_T_consistent_color(color)) { 
        is_have_tc = true; 
        break;
      }
    }

    if (is_have_tc) { 
      current.insert(edge.first); 
      processed.insert(edge.first.first); processed.insert(edge.first.second);
      processed.erase(Infty);
    } 
  } 

  std::set<std::set<edge_t> > answer;
  if (processed.size() == vertex_set.size()) { 
    answer.insert(current);
  } else { 
    std::function<void()> findLambda = [&] () -> void { 
      for (auto const &e : weight_edges) { 
        if (processed.count(e.first.first) == 0 && processed.count(e.first.second) == 0) { 
          processed.insert(e.first.first); processed.insert(e.first.second);
          processed.erase(Infty);
      	  current.insert(e.first);
          findLambda();
      	  current.erase(e.first);
      	  processed.erase(e.first.first);
      	  processed.erase(e.first.second);
        }
      }
      if (processed.size() == vertex_set.size()) { 
        answer.insert(current);
      } 
    }; 
 
    for (auto const &edge : weight_edges) { 
      if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) { 	
        current.insert(edge.first); 
        processed.insert(edge.first.first); processed.insert(edge.first.second);
        processed.erase(Infty);
        findLambda();  
        current.erase(edge.first);
        processed.erase(edge.first.first);
        processed.erase(edge.first.second);
      }
    }
  } 

  std::multimap<size_t, edge_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (auto const & match: answer) { 
    size_t weight = 0; 
    for (auto const & edge: match) { 
      weight += weight_edges.find(edge)->second.first;
    } 

    if (weight < min) { 
      min = weight; 
      minimal_matching.clear();
      for (auto const &edge: match) { 
        minimal_matching.insert(std::make_pair(weight_edges.find(edge)->second.first, edge));
      } 
    }

  } 
 
  /*
  TRACE("Find minimal matching");
  for (auto const & edge : minimal_matching) {
    TRACE("(" << edge.second.first  << ", " << edge.second.second << ") " << edge.first); 
  }
  */

  return minimal_matching;
} 

} 

#endif
