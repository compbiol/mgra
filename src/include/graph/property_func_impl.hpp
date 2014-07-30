#ifndef PROPERTY_FUNC_IMPL_HPP
#define PROPERTY_FUNC_IMPL_HPP

/*
 * Implementation of the function about property for vertex
 */

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_simple_vertex(vertex_t const & v) const {
  return (!is_duplication_vertex(v) && !is_indel_vertex(v) && (this->degree_vertex(v) == 2)); 
}  

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_have_self_loop(vertex_t const & v) const {
  return !(this->get_edge_multicolor(v, v).empty()); 
} 
 
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_indel_vertex(vertex_t const & v) const {
  if (is_duplication_vertex(v)) {
    return false; 
  } 
  mcolor_t const & union_color = this->get_adjacent_multiedges(v).union_multicolors();
  return (union_color.is_one_to_one_match() && (union_color != multicolors.get_complete_color()));
}

template<class mcolor_t>  
bool BreakpointGraph<mcolor_t>::is_duplication_vertex(vertex_t const & v) const {  
  if (is_have_self_loop(v)) {
    return true;
  } 
  mcolor_t const & union_color = this->get_adjacent_multiedges(v).union_multicolors();
  return !(union_color.is_one_to_one_match());
} 

/*
 * Implementation of the function about property for edges
 */
template<class mcolor_t>
size_t BreakpointGraph<mcolor_t>::calculate_cost(vertex_t const & u, vertex_t const & v) const {
  if (u == Infty) {
    mularcs_t mularcs_v = this->get_adjacent_multiedges_with_info(v, false);
    mularcs_v.erase(u);
    return mularcs_v.size();
  }

  if (v == Infty) {
    mularcs_t mularcs_u = this->get_adjacent_multiedges_with_info(u, false);
    mularcs_u.erase(v);
    return mularcs_u.size();
  }

  mularcs_t mularcs_v = this->get_adjacent_multiedges_with_info(v, false);
  mularcs_v.erase(u);
  mularcs_t mularcs_u = this->get_adjacent_multiedges_with_info(u, false);
  mularcs_u.erase(v);
    
  typedef std::tuple<vertex_t, mcolor_t, size_t> ind_acr_t; 
  utility::equivalence<ind_acr_t> equiv; 
  
  for (auto const & arc_u : mularcs_u) { 
    for (auto const & arc_v : mularcs_v) { 
      mcolor_t color(arc_u.second, arc_v.second, mcolor_t::Intersection);
      if (color.size() > 0) { 
        equiv.addrel(ind_acr_t(arc_u.first, arc_u.second, 0), ind_acr_t(arc_v.first, arc_v.second, 1));
      } 
    } 
  }
  
  equiv.update();
  std::map<ind_acr_t, std::set<ind_acr_t> > const & classes = equiv.template get_eclasses<std::set<ind_acr_t> >(); 

  size_t count_U = 0; 
  for (auto const & color_set : classes) { 
    count_U += (color_set.second.size() - 1);
  }
  return count_U;
}

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const {
  if (this->postponed_deletions.defined(x, y)) {
    return false;
  } 

  if (!multicolors.is_vec_T_consistent_color(color)) {
    return false;
  } 

  mularcs_t mularcs_x;
  if (x != Infty) {
    mularcs_x = this->get_adjacent_multiedges_with_info(x);
    mularcs_x.erase(y);
  } 

  mularcs_t mularcs_y;
  if (y != Infty) {
    mularcs_y = this->get_adjacent_multiedges_with_info(y);
    mularcs_y.erase(x);
  } 

  bool mobile = false;

  for(auto arc_x = mularcs_x.cbegin(); (arc_x != mularcs_x.cend()) && !mobile; ++arc_x) { 
    mobile = canformQ(arc_x->first, color);
  } 

  for(auto arc_y = mularcs_y.cbegin(); (arc_y != mularcs_y.cend()) && !mobile; ++arc_y) { 
    mobile = canformQ(arc_y->first, color);
  } 

  return mobile;    
}

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_mobility_edge(vertex_t const & x, vertex_t const & y) const {
  if (this->postponed_deletions.defined(x, y)) {
    return false;
  } 
  
  if ((x == Infty || y == Infty) && is_mobile_irregular_edge) { 
    return true;
  }

  bool non_mobile = false; 
  auto const & central_colors = this->get_edge_multicolor_with_info(x, y);
  for (auto color = central_colors.cbegin(); (color != central_colors.cend()) && !non_mobile; ++color) { 
    if (!multicolors.is_T_consistent_color(*color)) { 
      non_mobile = true;
    } else if (!this->is_vec_T_consistent_color(*color)) { 
      non_mobile = true;
    } 
  }

  if (non_mobile) { 
    return false; 
  }

  mularcs_t mularcs_x;
  if (x != Infty) { 
    mularcs_x = this->get_adjacent_multiedges_with_info(x);
    mularcs_x.erase(y);
  } 

  mularcs_t mularcs_y;
  if (y != Infty) {
    mularcs_y = this->get_adjacent_multiedges_with_info(y); 
    mularcs_y.erase(x); 
  } 

  bool mobile = false;  
  for (auto color = central_colors.cbegin(); (color != central_colors.cend()) && !mobile; ++color) { 
    for(auto arc_x = mularcs_x.cbegin(); (arc_x != mularcs_x.cend()) && !mobile; ++arc_x) { 
      mobile = this->canformQ(arc_x->first, *color);
    }

    for(auto arc_y = mularcs_y.cbegin(); (arc_y != mularcs_y.cend()) && !mobile; ++arc_y) { 
      mobile = this->canformQ(arc_y->first, *color);
    } 
  }

  return mobile; 
} 

#endif