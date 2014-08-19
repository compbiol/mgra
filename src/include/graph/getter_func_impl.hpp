#ifndef GETTER_FUNC_IMPL_HPP
#define GETTER_FUNC_IMPL_HPP

/*
 * Implementation of the function to get all incident edges from vertex v. 
 */
template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_all_adjacent_multiedges(vertex_t const & u) const {  
  mularcs_t current = graph.get_all_adjacent_multiedges<mularcs_t>(u);
  
#ifdef PSEUDO_EDGE
  mularcs_t split; 
  for (auto const & arc : current) { 
    auto colors = split_edge_on_pseudo_edge(u, arc.second, arc.first); 
    for (auto const & color : colors) {  
      split.insert(arc.first, color); 
    }
  }
  return split;
#else
  return current;
#endif 

}  

template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_all_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge) const { 
  mularcs_t output = get_all_adjacent_multiedges(u);
  
  if (this->number_of_splits != 1) { 
    mularcs_t split; 
    for(auto const & arc : output) {
      if ((!with_bad_edge || (with_bad_edge && !this->postponed_deletions.defined(u, arc.first))) 
        && !multicolors.is_T_consistent_color(arc.second) && arc.second.size() < graph.count_local_graphs()) {
        auto const & colors = multicolors.split_color_on_tc_color(arc.second, number_of_splits);
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

#ifdef PSEUDO_EDGE
template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_real_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge) const { 
  mularcs_t const & current = get_all_adjacent_multiedges_with_info(u, with_bad_edge);

  mularcs_t output; 

  for (auto const & arc : current) { 
    if (!is_pseudo_edge(u, arc.second)) { 
      output.insert(arc.first, arc.second);
    } 
  } 
  
  return output;
} 

template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_clones_adjacent_multiedges_with_info(vertex_t const & u, mcolor_t const & color, bool with_bad_edge) const { 
  mularcs_t const & mularcs_u = get_all_adjacent_multiedges_with_info(u, with_bad_edge);

  assert(mother_verteces.count(u) != 0);

  auto const & clones = this->mother_verteces.find(u)->second;
  bool flag = false;
  size_t current_clone = 0;

  for (auto clone_ind = clones.crbegin(); clone_ind != clones.crend() && !flag; ++clone_ind) { 
    mcolor_t const & clone_color = history.get_clone(*clone_ind).get_mcolor();
    mcolor_t inter_color(clone_color, color, mcolor_t::Intersection);
    if (inter_color.size() > 0 && inter_color.size() < clone_color.size()) { 
      current_clone = *clone_ind;
      flag = true;
    }
  }

  mularcs_t output; 

  if (flag) {  
    for (auto const & arc : mularcs_u) {
      mcolor_t const & clone_color = history.get_clone(current_clone).get_mcolor();
      mcolor_t inter_color(arc.second, clone_color, mcolor_t::Intersection);
      if (inter_color.size() > 0 && inter_color.size() < clone_color.size()) { 
        auto edges = history.get_clone(current_clone).get_end_edges();
        for (auto const & clone_arc : edges) { 
          if (arc.second.includes(clone_arc.second)) { 
            output.insert(arc.first, arc.second);
            break;
          }    
        } 
      } 
    } 
  } 

  return output;
}
#endif

template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_all_adjacent_vtc_multiedges(vertex_t const & u) const { 
  mularcs_t current = graph.get_all_adjacent_multiedges<mularcs_t>(u);
  
  mularcs_t split;
  for(auto const & arc : current) {
    if (!multicolors.is_vec_T_consistent_color(arc.second) && arc.second.size() < graph.count_local_graphs()) {
      auto const & colors = multicolors.split_color_on_vtc_color(arc.second);
      for(auto const & color : colors) {  
        split.insert(arc.first, color); 
      }
    } else { 
      split.insert(arc.first, arc.second); 
    }
  } 

  return split; 
}

template<class mcolor_t>
structure::Mularcs<mcolor_t> BreakpointGraph<mcolor_t>::get_all_adjacent_tc_multiedges(vertex_t const & u) const { 
  mularcs_t current = get_all_adjacent_multiedges(u);
  
  mularcs_t split;

  for(auto const & arc : current) {
    if (!multicolors.is_T_consistent_color(arc.second) && arc.second.size() < graph.count_local_graphs()) {
      auto const & colors = multicolors.split_color_on_tc_color(arc.second, graph.count_local_graphs() + 1);
      for(auto const & color : colors) {  
        split.insert(arc.first, color); 
      }
    } else { 
      split.insert(arc.first, arc.second); 
    }
  } 

  return split; 
}

/*
 * Implementation of the function to get multiocolor for edge (u, v). 
 */
template<class mcolor_t>
mcolor_t BreakpointGraph<mcolor_t>::get_all_multicolor_edge(vertex_t const & u, vertex_t const & v) const { 
  return graph.get_all_multicolor_edge<mcolor_t>(u, v);
}

#ifdef PSEUDO_EDGE
template<class mcolor_t>
mcolor_t BreakpointGraph<mcolor_t>::get_real_multicolor_edge(vertex_t const & x, vertex_t const & y) const { 
  mcolor_t current = graph.get_all_multicolor_edge<mcolor_t>(x, y);
  std::set<mcolor_t> output = split_edge_on_pseudo_edge(x, current, y); 
  
  mcolor_t result;

  for (auto const & color : output) { 
    if (!is_pseudo_edge(x, color, y)) { 
      result = mcolor_t(color, result, mcolor_t::Union);
    } 
  } 

  return result;
}
#endif 

template<class mcolor_t>
std::set<mcolor_t> BreakpointGraph<mcolor_t>::get_all_multicolor_edge_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge) const { 
  assert(u != Infty || v != Infty);

  mcolor_t current = get_all_multicolor_edge(u, v);
#ifdef PSEUDO_EDGE
  std::set<mcolor_t> output = split_edge_on_pseudo_edge(u, current, v); 
#else
  std::set<mcolor_t> output({current});
#endif

  if (this->number_of_splits != 1) { 
    std::set<mcolor_t> results;
#ifdef PSEUDO_EDGE
    for (auto const & color : output) { 
      if ((!with_bad_edge || (with_bad_edge && !postponed_deletions.defined(u, v))) 
        && !multicolors.is_T_consistent_color(color) && color.size() < graph.count_local_graphs()) {
        auto temp = multicolors.split_color_on_tc_color(color, number_of_splits);
        results.insert(temp.begin(), temp.end());
      } else { 
        results.insert(color);
      }
    }
#else 
    if ((!with_bad_edge || (with_bad_edge && !this->postponed_deletions.defined(u, v))) 
        && !multicolors.is_T_consistent_color(current) && current.size() < graph.count_local_graphs()) {
      auto temp = multicolors.split_color_on_tc_color(current, number_of_splits);
      results.insert(temp.begin(), temp.end());
    } else { 
      results.insert(current);
    }
#endif
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

  for(vertex_t const & x : graph) {
    if (!not_drop_complete_edge) {
      connected_components.addrel(x, x);
    } 
    
    mularcs_t const & mularcs = get_all_adjacent_multiedges(x); 
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

  for(vertex_t const & x : graph) {
    mularcs_t const & mularcs = get_all_adjacent_multiedges(x); 
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

/*
 * Implementation of the function to get multiocolor for edge (u, v). 
 */
#ifdef PSEUDO_EDGE
template<class mcolor_t>
std::pair<vertex_t, vertex_t> BreakpointGraph<mcolor_t>::get_real_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const {

  auto const & change_lambda = [&](vertex_t const & a) -> std::pair<bool, vertex_t> { 
    std::pair<bool, vertex_t> result(false, a);
    
    if (this->mother_verteces.count(a) != 0) {
      auto const & clones = this->mother_verteces.find(a)->second;
      for (auto clone_ind = clones.crbegin(); clone_ind != clones.crend() && !result.first; ++clone_ind) { 
        auto clone = history.get_clone(*clone_ind);
        mcolor_t const & clone_color = clone.get_mcolor();
        mcolor_t inter_color(clone_color, color, mcolor_t::Intersection);
        if (inter_color.size() > 0 && inter_color.size() < clone_color.size() && inter_color.size() != color.size()) { 
          result = std::make_pair(true, clone.get_father());  
        }
      }
    }

    return result; 
  };         

  std::pair<bool, vertex_t> left(true, x);
  while (left.first) { 
    left = change_lambda(left.second);
  }

  std::pair<bool, vertex_t> right(true, y);
  while (right.first) { 
    right = change_lambda(right.second);
  }

  return std::make_pair(left.second, right.second);
}
#endif

#endif
