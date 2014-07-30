#ifndef GETTER_FUNC_IMPL_HPP
#define GETTER_FUNC_IMPL_HPP

/*
 * Implementation of the function to get all incident edges from vertex v. 
 */
template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_adjacent_multiedges(vertex_t const & u) const { 
  assert(u != Infty);
  
  mularcs_t output;

  for (size_t i = 0; i < this->m_local_graphs.size(); ++i) {
    auto iters = this->m_local_graphs[i].equal_range(u);
    for (auto it = iters.first; it != iters.second; ++it) { 
      output.insert(it->second, i); 
    }
  }

  return output;
  /*mularcs_t split; 
  for (auto const & arc : output) { 
    auto colors = split_edge_on_pseudo_edge(u, arc.second, arc.first); 
    for (auto const & color : colors) {  
      split.insert(arc.first, color); 
    }
  }
  return split;*/
} 

template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge) const { 
  assert(u != Infty);

  mularcs_t output = this->get_adjacent_multiedges(u);
  
  if (this->number_of_splits != 1) { 
    mularcs_t split; 
    for(auto const & arc : output) {
      if ((!with_bad_edge || (with_bad_edge && !this->postponed_deletions.defined(u, arc.first))) 
        && !multicolors.is_T_consistent_color(arc.second) && arc.second.size() < this->count_local_graphs()) {
        auto const & colors = multicolors.split_color(arc.second, number_of_splits);
        for(auto const & color : colors) {  
          split.insert(arc.first, color); 
        }
      } else { 
        split.insert(arc.first, arc.second); 
      }
    }
    return split; 
  }

  return output;
} 

/*
 * Implementation of the function to get multiocolor for edge (u, v). 
 */
template<class mcolor_t>
mcolor_t BreakpointGraph<mcolor_t>::get_edge_multicolor(vertex_t const & u, vertex_t const & v) const { 
  assert(u != Infty || v != Infty);

  mcolor_type result;

  for (size_t i = 0; i < this->m_local_graphs.size(); ++i) {
    auto iters = this->m_local_graphs[i].equal_range(u);
    for (auto it = iters.first; it != iters.second; ++it) { 
      if (it->second == v) { 
        result.insert(i);
      }
    } 
  } 

  return result;
}

template<class mcolor_t>
std::set<mcolor_t> BreakpointGraph<mcolor_t>::get_edge_multicolor_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge) const { 
  assert(u != Infty || v != Infty);

  mcolor_t current = get_edge_multicolor(u, v);
  std::set<mcolor_t> output({current});// = split_edge_on_pseudo_edge(u, current, v);

  if (this->number_of_splits != 1) { 
    std::set<mcolor_t> results;
    if ((!with_bad_edge || (with_bad_edge && !this->postponed_deletions.defined(u, v))) 
        && !multicolors.is_T_consistent_color(current) && current.size() < this->count_local_graphs()) {
      auto temp = multicolors.split_color(current, number_of_splits);
      results.insert(temp.begin(), temp.end());
    } else { 
      results.insert(current);
    }
    //for (auto const & color : output) { 
      /*if ((!with_bad_edge || (with_bad_edge && !postponed_deletions.defined(u, v))) 
        && !multicolors.is_T_consistent_color(color) && color.size() < count_local_graphs()) {
        auto temp = split_color(color, number_of_splits);
        results.insert(temp.begin(), temp.end());
      } else { 
        results.insert(color);
      }*/
    //} 
    return results;
  } 

  return output;
}


/*
 * Implementation of the function for split graph on components
 */
template<class mcolor_t>
utility::equivalence<vertex_t> BreakpointGraph<mcolor_t>::split_on_components(bool not_drop_complete_edge) const { 
  utility::equivalence<vertex_t> connected_components; // connected components

  for(vertex_t const & x : *this) {
    if (!not_drop_complete_edge) {
      connected_components.addrel(x, x);
    } 
    
    mularcs_t const & mularcs = get_adjacent_multiedges(x); 
    if (not_drop_complete_edge && mularcs.size() == 1 && mularcs.union_multicolors() == multicolors.get_complete_color()) { 
      continue; // ignore complete multiedges
    } 

    for(auto const & arc : mularcs) {    
      if (arc.first != Infty) { 
        connected_components.addrel(x, arc.first);
      } 
    }
  }
        
  connected_components.update();
   
  return connected_components;
}

template<class mcolor_t>
utility::equivalence<vertex_t> BreakpointGraph<mcolor_t>::split_on_components_with_color(mcolor_t const & color) const { 
  utility::equivalence<vertex_t> connected_components; // connected components

  for(vertex_t const & x : *this) {
    mularcs_t const & mularcs = get_adjacent_multiedges(x); 
    if (mularcs.size() == 1) { 
      continue;
    }

    for(auto const & arc : mularcs) {    
      if (arc.first != Infty && arc.second != color) { 
        connected_components.addrel(x, arc.first);
      }
    }
  } 

  connected_components.update();

  return connected_components; 
}


 #endif