#ifndef GETTER_FUNC_IMPL_HPP
#define GETTER_FUNC_IMPL_HPP

/*
 * Implementation of the function to get all incident edges from vertex u. 
 */
template<class mcolor_t>
structure::Mularcs<mcolor_t> GraphPack<mcolor_t>::get_all_adjacent_multiedges(vertex_t const & u) const {  
  mularcs_t current = graph.get_all_adjacent_multiedges<mularcs_t>(u);  
  return current;
}  

template<class mcolor_t>
structure::Mularcs<mcolor_t> GraphPack<mcolor_t>::get_all_adjacent_multiedges_with_info(vertex_t const & u, bool with_bad_edge) const { 
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

template<class mcolor_t>
structure::Mularcs<mcolor_t> GraphPack<mcolor_t>::get_all_adjacent_tc_multiedges(vertex_t const & u) const { 
  mularcs_t current = graph.get_all_adjacent_multiedges<mularcs_t>(u);
  
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

template<class mcolor_t>
structure::Mularcs<mcolor_t> GraphPack<mcolor_t>::get_all_adjacent_vtc_multiedges(vertex_t const & u) const { 
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

/*
 * Implementation of the function to get multiocolor for edge (u, v). 
 */
template<class mcolor_t>
mcolor_t GraphPack<mcolor_t>::get_all_multicolor_edge(vertex_t const & u, vertex_t const & v) const { 
  return graph.get_all_multicolor_edge<mcolor_t>(u, v);
}

template<class mcolor_t>
std::set<mcolor_t> GraphPack<mcolor_t>::get_all_multicolor_edge_with_info(vertex_t const & u, vertex_t const & v, bool with_bad_edge) const { 
  assert(u != Infty || v != Infty);

  mcolor_t current = get_all_multicolor_edge(u, v);
  std::set<mcolor_t> output;
  output.insert(current);

  if (this->number_of_splits != 1) { 
    std::set<mcolor_t> results;

    if ((!with_bad_edge || (with_bad_edge && !this->postponed_deletions.defined(u, v))) 
        && !multicolors.is_T_consistent_color(current) && current.size() < graph.count_local_graphs()) {
      auto temp = multicolors.split_color_on_tc_color(current, number_of_splits);
      results.insert(temp.begin(), temp.end());
    } else { 
      results.insert(current);
    }

    return results;
  } 

  return output;
}

/*
 * Implementation of the function for split graph on components
 */
template<class mcolor_t>
utility::equivalence<vertex_t> GraphPack<mcolor_t>::split_on_components(bool not_drop_complete_edge) const { 
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
utility::equivalence<vertex_t> GraphPack<mcolor_t>::split_on_components_with_color(mcolor_t const & color) const { 
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

#endif
