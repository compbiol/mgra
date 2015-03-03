#ifndef APPLY_FUNCTION_IMPL_HPP
#define APPLY_FUNCTION_IMPL_HPP

/*
 * Implementation of the function for operation on graph, which we do. 
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::apply(twobreak_t const & twobreak, bool record) { 
  assert(multicolors.is_vec_T_consistent_color(twobreak.get_mcolor()));

  if (record) {
    history.save_twobreak(twobreak);
  }
 
  for (auto const & color : twobreak) {
    for (size_t i = 0; i < 2; ++i) {
      if (twobreak.get_arc(i).first != Infty || twobreak.get_arc(i).second != Infty) {
        graph.erase_edge(color.first, twobreak.get_arc(i).first, twobreak.get_arc(i).second);
      } 
    }

    if (twobreak.get_vertex(0) != Infty || twobreak.get_vertex(2) != Infty) {
      graph.add_edge(color.first, twobreak.get_vertex(0), twobreak.get_vertex(2));
    }

    if (twobreak.get_vertex(1) != Infty || twobreak.get_vertex(3) != Infty) {
      graph.add_edge(color.first, twobreak.get_vertex(1), twobreak.get_vertex(3));
    }
  }
  
  check_changed_vertex(twobreak.get_vertex(0)); 
  check_changed_vertex(twobreak.get_vertex(1)); 
  check_changed_vertex(twobreak.get_vertex(2)); 
  check_changed_vertex(twobreak.get_vertex(3));
} 

template<class mcolor_t>
void GraphPack<mcolor_t>::apply(clone_t const & clone, bool record) { 
  size_t number; 

  if (record) { 
    number = history.save_clone(clone);
  } 

  auto const & central_edge = clone.get_central_arc(); 
  auto const & mother_edge = clone.get_mother_edge();

  assert(multicolors.is_vec_T_consistent_color(mother_edge.second));

  if (clone.is_have_pseudo_vertex()) {
    graph.add_vertex(mother_edge.first);
    for (auto const & color : mother_edge.second) {
      graph.erase_edge(color.first, central_edge.second, Infty);  
      graph.add_edge(color.first, central_edge.second, mother_edge.first);
    }  

    mcolor_t const & compl_color = this->multicolors.get_complement_color(mother_edge.second);
    for (auto const & color : compl_color) {
      graph.add_edge(color.first, mother_edge.first, Infty);
    }
  }

  for (auto const & color : mother_edge.second) {
    graph.erase_edge(color.first, central_edge.second, mother_edge.first);  
    graph.add_edge(color.first, central_edge.first, central_edge.second);
  } 

  mularcs_t const & fathers = clone.get_end_edges();
  for (auto const & arc: fathers) { 
    for (auto const & color : arc.second) {
      graph.erase_edge(color.first, arc.first, central_edge.first);  
      graph.add_edge(color.first, arc.first, mother_edge.first);
    } 
  } 

  this->mother_verteces[mother_edge.first].push_back(number);

  check_changed_vertex(mother_edge.first);
  check_changed_vertex(central_edge.first);
  check_changed_vertex(central_edge.second);
  for (auto const & arc: fathers) { 
    check_changed_vertex(arc.first);
  } 
} 

template<class mcolor_t>
void GraphPack<mcolor_t>::apply(insdel_t const & insdel, bool record) { 
  if (record) { 
    ;
  }

  for (auto const &color : insdel) {
    graph.add_edge(color.first, insdel.get_edge().first, insdel.get_edge().second); 
  }

  if (insdel.is_insertion()) {
    prosthetic_chromosomes.insert(insdel.get_edge().first, insdel.get_edge().second);
  } else { 
    prosthetic_chromosomes.insert(insdel.get_edge().first, insdel.get_edge().second);
    postponed_deletions.insert(insdel.get_edge().first, insdel.get_edge().second);
  }   
}

template<class mcolor_t>
void GraphPack<mcolor_t>::apply(tandem_duplication_t const & duplication, bool record) {
  if (record) {
    history.save_tandem_duplication(duplication);
  }
 
  if (duplication.is_deletion_oper()) {
    for (auto it = duplication.cbegin_edges(); it != (--duplication.cend_edges()); ++it) {
      for (auto const &color : duplication) {  
        graph.erase_edge(color.first, it->first, it->second);
      }
    } 

    if (duplication.is_reverse_tandem_duplication()) {
      for (auto const &color : duplication) {
        graph.add_edge(color.first, (--duplication.cend_edges())->first, (--duplication.cend_edges())->second);
      }
    } else {
      for (auto const &color : duplication) {
        graph.erase_edge(color.first, (--duplication.cend_edges())->first, (--duplication.cend_edges())->second);
      }
    }
  } else {
    for (auto it = duplication.cbegin_edges(); it != (--duplication.cend_edges()); ++it) {
      for (auto const &color : duplication) {
        graph.add_edge(color.first, it->first, it->second);
      }
    } 

    if (duplication.is_reverse_tandem_duplication()) {
      for (auto const &color : duplication) {
        graph.erase_edge(color.first, (--duplication.cend_edges())->first, (--duplication.cend_edges())->second);
      }
    } else {
      for (auto const &color : duplication) {
        graph.add_edge(color.first, (--duplication.cend_edges())->first, (--duplication.cend_edges())->second);
      }
    }
  } 
} 

#endif