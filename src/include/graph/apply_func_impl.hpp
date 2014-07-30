#ifndef APPLY_FUNCTION_IMPL_HPP
#define APPLY_FUNCTION_IMPL_HPP

/*
 * Implementation of the function for operation on graph, which we do. 
 */

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::apply(twobreak_t const & twobreak, bool record) { 
  assert(multicolors.is_vec_T_consistent_color(twobreak.get_mcolor()));

  if (record) {
    history.save_twobreak(twobreak);
  }
 
  for (auto const & color : twobreak) {
    for (size_t i = 0; i < 2; ++i) {
      if (twobreak.get_arc(i).first != Infty || twobreak.get_arc(i).second != Infty) {
        this->erase_edge(color.first, twobreak.get_arc(i).first, twobreak.get_arc(i).second);
      } 
    }

    if (twobreak.get_vertex(0) != Infty || twobreak.get_vertex(2) != Infty) {
      this->add_edge(color.first, twobreak.get_vertex(0), twobreak.get_vertex(2));
    }

    if (twobreak.get_vertex(1) != Infty || twobreak.get_vertex(3) != Infty) {
      this->add_edge(color.first, twobreak.get_vertex(1), twobreak.get_vertex(3));
    }
  }

  
  check_changed_vertex(twobreak.get_vertex(0)); 
  check_changed_vertex(twobreak.get_vertex(1)); 
  check_changed_vertex(twobreak.get_vertex(2)); 
  check_changed_vertex(twobreak.get_vertex(3));
} 

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::apply(clone_t const & clone, bool record) { 
  size_t number; 

  if (record) { 
    number = history.save_clone(clone);
  } 

  auto const & central_edge = clone.get_central_arc(); 
  auto const & mother_edge = clone.get_mother_edge();

  assert(multicolors.is_vec_T_consistent_color(mother_edge.second));

  if (clone.is_have_pseudo_vertex()) {
    this->add_vertex(mother_edge.first);
    for (auto const & color : mother_edge.second) {
      this->erase_edge(color.first, central_edge.second, Infty);  
      this->add_edge(color.first, central_edge.second, mother_edge.first);
    }  

    mcolor_t const & compl_color = this->get_complement_color(mother_edge.second);
    for (auto const & color : compl_color) {
      this->add_edge(color.first, mother_edge.first, Infty);
    }
  }

  for (auto const & color : mother_edge.second) {
    this->erase_edge(color.first, central_edge.second, mother_edge.first);  
    this->add_edge(color.first, central_edge.first, central_edge.second);
  } 

  mularcs_t const & fathers = clone.get_end_edges();
  for (auto const & arc: fathers) { 
    for (auto const & color : arc.second) {
      this->erase_edge(color.first, arc.first, central_edge.first);  
      this->add_edge(color.first, arc.first, mother_edge.first);
    } 
  } 
 
  this->pseudo_edges[mother_edge.first].insert(mother_edge.second);
  if (clone.is_have_pseudo_vertex() && record) {
    this->pseudo_infinity_verteces[mother_edge.first].insert(std::make_pair(mother_edge.second, number));
  } 

  check_changed_vertex(mother_edge.first);
  check_changed_vertex(central_edge.first);
  check_changed_vertex(central_edge.second);
  for (auto const & arc: fathers) { 
    check_changed_vertex(arc.first);
  } 
} 

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::apply(insertion_t const & insdel, bool record) { 
  if (record) { 
    ;
  }

  for (auto const &color : insdel) {
    this->add_edge(color.first, insdel.get_edge().first, insdel.get_edge().second); 
  }

  if (insdel.is_insertion()) {
    this->insertions.insert(insdel.get_edge().first, insdel.get_edge().second);
  } else { 
    this->postponed_deletions.insert(insdel.get_edge().first, insdel.get_edge().second);
  }   
}

template<class mcolor_t>
void BreakpointGraph<mcolor_t>::apply(tandem_duplication_t const & dupl, bool record) {
  if (record) {
    history.save_tandem_duplication(dupl);
  }
 
  if (dupl.is_deletion_oper()) {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (auto const &color : dupl) {  
        this->erase_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (auto const &color : dupl) {
        this->add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (auto const &color : dupl) {
        this->erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } else {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (auto const &color : dupl) {
        this->add_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (auto const &color : dupl) {
        this->erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (auto const &color : dupl) {
        this->add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } 
} 

#endif