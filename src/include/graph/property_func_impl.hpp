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
  return !(this->get_all_multicolor_edge(v, v).empty()); 
} 
 
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_indel_vertex(vertex_t const & v) const {
  if (v == Infty) { 
    return false;
  }
  if (is_duplication_vertex(v)) {
    return false; 
  } 
  mcolor_t const & union_color = this->get_all_adjacent_multiedges(v).union_multicolors();
  return (union_color.is_one_to_one_match() && (union_color != multicolors.get_complete_color()));
}

template<class mcolor_t>  
bool BreakpointGraph<mcolor_t>::is_duplication_vertex(vertex_t const & v) const {  
  if (v == Infty) { 
    return false;
  }
  if (is_have_self_loop(v)) {
    return true;
  } 
  mcolor_t const & union_color = this->get_all_adjacent_multiedges(v).union_multicolors();
  return !(union_color.is_one_to_one_match());
} 

/*
 * Implementation of the function about property for edges
 */
template<class mcolor_t>
size_t BreakpointGraph<mcolor_t>::mobility_score(graph_edge_t const & viewed, mcolor_t const & color, graph_edge_t const & removed) const {
  auto calculate_lambda = [&] (vertex_t const & target, vertex_t const & other) -> size_t {
    mularcs_t mularcs; 
#ifdef PSEUDO_EDGE
    if (is_pseudo_edge(target, color)) { 
      mularcs = get_clones_adjacent_multiedges_with_info(target, color);
    } else if (target != Infty) { 
      mularcs = get_real_adjacent_multiedges_with_info(target);
    }
#else
    if (target != Infty) { 
      mularcs = get_all_adjacent_multiedges_with_info(target); 
    } 
#endif

    mularcs.erase(other);
    mularcs.erase(removed.first); 
    mularcs.erase(removed.second); 
    
    size_t number_variant = 0; 
    for (auto const & arc : mularcs) { 
      if (this->canformQ(arc, color)) {
        ++number_variant;
      } 
    }
    return number_variant;
  };

  return (calculate_lambda(viewed.first, viewed.second) + calculate_lambda(viewed.second, viewed.first)); 
}

#ifdef PSEUDO_EDGE
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_result_mobility_edge(graph_edge_t const & viewed, mcolor_t const & new_color, graph_edge_t const & removed) const { 
  graph_edge_t result_viewed(viewed.first, removed.first);
  graph_edge_t result_another(viewed.second, removed.second);

  auto const is_have_multiedge_lambda = [&](vertex_t const & x, vertex_t const & y) -> bool {
    mcolor_t union_color(get_all_multicolor_edge(x, y), new_color, mcolor_t::Union);   
    std::set<mcolor_t> unions = split_edge_on_pseudo_edge(x, union_color, y); 
  
    std::set<mcolor_t> results; 

    for (auto const & color : unions) { 
      if (number_of_splits != 1 && !postponed_deletions.defined(x, y) 
          && !multicolors.is_T_consistent_color(color) && color.size() < graph.count_local_graphs()) {
        auto temp = multicolors.split_color_on_tc_color(color, number_of_splits);
        results.insert(temp.begin(), temp.end());
      } else { 
        results.insert(color);
      }
    }

    return (results.count(new_color) != 0);
  };  
  
  if (is_have_multiedge_lambda(result_viewed.first, result_viewed.second)) { 
    mularcs_t mularcs_first; 
    if (this->is_pseudo_edge(result_viewed.first, new_color)) {
      mularcs_first = get_clones_adjacent_multiedges_with_info(result_viewed.first, new_color); 
    } else if (result_viewed.first != Infty) { 
      mularcs_first = get_real_adjacent_multiedges_with_info(result_viewed.first);
    }
    mularcs_first.erase(result_viewed.second);

    bool is_avaliable = false; 
    for (auto arc = mularcs_first.cbegin(); arc != mularcs_first.cend() && !is_avaliable; ++arc) {
      if ((arc->first == result_another.second || arc->first == result_another.first) 
          && arc->second != new_color) {
        is_avaliable = is_have_multiedge_lambda(result_another.first, result_another.second);
      }
    }

    mularcs_t mularcs_second;
    if (this->is_pseudo_edge(result_viewed.second, new_color)) {
      mularcs_second = get_clones_adjacent_multiedges_with_info(result_viewed.second, new_color); 
    } else if (result_viewed.second != Infty) { 
      mularcs_second = get_real_adjacent_multiedges_with_info(result_viewed.second);
    }
    mularcs_second.erase(result_viewed.first);

    for (auto arc = mularcs_first.cbegin(); arc != mularcs_first.cend() && !is_avaliable; ++arc) {
      if ((arc->first == result_another.second || arc->first == result_another.first) 
          && arc->second != new_color) {
        is_avaliable = is_have_multiedge_lambda(result_another.first, result_another.second);
      }
    }

    return is_avaliable;
  } 

  return false;
}
#endif


template<class mcolor_t>
std::pair<size_t, size_t> BreakpointGraph<mcolor_t>::is_decrease_verteces_score(twobreak_t const & twobreak) const { 
  auto calc_score_lambda = [&] (mularcs_t & mularcs) -> size_t { 
    vertex_t target_vertex; 
    size_t score_max = 0; 
    std::unordered_map<vertex_t, size_t> how_many; 
    for (auto const & arc : mularcs) { 
      if (!multicolors.is_vec_T_consistent_color(arc.second)) { 
        assert ((multicolors.is_T_consistent_color(arc.second)) || (arc.second == multicolors.get_complete_color()));
        target_vertex = arc.first; 
        score_max = graph.count_local_graphs(); 
        break;
      } else { 
        auto iter = how_many.find(arc.first);
        if (iter != how_many.end()) { 
          iter->second = 1 + iter->second;
        } else { 
          how_many.insert(std::make_pair(arc.first, 1));
        } 

        if (how_many.find(arc.first)->second > score_max) {
          score_max = how_many.find(arc.first)->second;
          vertex_t target_vertex = arc.first; 
        }
      }
    }   

    mularcs.erase(target_vertex); 

    return mularcs.size(); 
  };

  auto calc_vertex_score_lambda = [&] (vertex_t const & vertex) -> size_t { 
    mularcs_t mularcs; 
    if (vertex != Infty) {
      mularcs = this->get_all_adjacent_tc_multiedges(vertex);
    } 
    return calc_score_lambda(mularcs);
  };

  auto calc_future_score_lambda = [&] (vertex_t const & vertex, mcolor_t const & color, vertex_t const & removed, vertex_t const & inserted) -> size_t { 
    if (vertex != Infty) {
      mularcs_t mularcs = this->get_all_adjacent_tc_multiedges(vertex);
      /*if (vertex == "380t" || vertex == "793t") { 
        std::cerr << "Init " << vertex << std::endl;
        std::cerr << genome_match::mcolor_to_name(color) << std::endl;
        for (auto const & arc : mularcs) {
          std::cerr << arc.first << " " << genome_match::mcolor_to_name(arc.second) << std::endl;
        }
      } */

      mularcs.erase(removed, color); 

      /*if (vertex == "380t" || vertex == "793t") { 
        std::cerr << "after erase " << vertex << std::endl;
        for (auto const & arc : mularcs) {
          std::cerr << arc.first << " " << genome_match::mcolor_to_name(arc.second) << std::endl;
        }
      } */

      auto all_colors = mularcs.equal_range(inserted);
      mcolor_t union_color = color; 
      for (auto iter = all_colors.first; iter != all_colors.second; ++iter) { 
        union_color = mcolor_t(union_color, iter->second, mcolor_t::Union);
      }
      mularcs.erase(inserted); 

      std::set<mcolor_t> results = multicolors.split_color_on_tc_color(union_color, graph.count_local_graphs() + 1);          
      for (mcolor_t const & col : results) {
        mularcs.insert(inserted, col);
      }
     
      /*if (vertex == "380t" || vertex == "793t") { 
        std::cerr << "result for " << vertex << std::endl;
        for (auto const & arc : mularcs) {
          std::cerr << arc.first << " " << genome_match::mcolor_to_name(arc.second) << std::endl;
        }
      }*/
      return calc_score_lambda(mularcs);
    } else { 
      return 0; 
    }
  };

  size_t before_score = calc_vertex_score_lambda(twobreak.get_vertex(0)) + calc_vertex_score_lambda(twobreak.get_vertex(1))
        + calc_vertex_score_lambda(twobreak.get_vertex(2)) + calc_vertex_score_lambda(twobreak.get_vertex(3));

  /*if (twobreak.get_vertex(0) == "380t" && twobreak.get_vertex(1) == "604h" && twobreak.get_vertex(2) == "793t" && twobreak.get_vertex(3) == "223h") {
    std::cerr << "Start for inverse twobreak" << std::endl;    
    std::cerr << calc_vertex_score_lambda(twobreak.get_vertex(0)) << std::endl; 
    std::cerr << calc_vertex_score_lambda(twobreak.get_vertex(1)) << std::endl;
    std::cerr << calc_vertex_score_lambda(twobreak.get_vertex(2)) << std::endl;
    std::cerr << calc_vertex_score_lambda(twobreak.get_vertex(3)) << std::endl;
    std::cerr << "Futture " << std::endl;
    std::cerr << calc_future_score_lambda(twobreak.get_vertex(0), twobreak.get_mcolor(), twobreak.get_vertex(1), twobreak.get_vertex(2)) << std::endl; 
    std::cerr << calc_future_score_lambda(twobreak.get_vertex(1), twobreak.get_mcolor(), twobreak.get_vertex(0), twobreak.get_vertex(3)) << std::endl; 
    std::cerr << calc_future_score_lambda(twobreak.get_vertex(2), twobreak.get_mcolor(), twobreak.get_vertex(3), twobreak.get_vertex(0)) << std::endl; 
    std::cerr << calc_future_score_lambda(twobreak.get_vertex(3), twobreak.get_mcolor(), twobreak.get_vertex(2), twobreak.get_vertex(1)) << std::endl; 
  }*/

  size_t future_score = calc_future_score_lambda(twobreak.get_vertex(0), twobreak.get_mcolor(), 
                  twobreak.get_vertex(1), twobreak.get_vertex(2));
  future_score += calc_future_score_lambda(twobreak.get_vertex(1), twobreak.get_mcolor(), 
                  twobreak.get_vertex(0), twobreak.get_vertex(3));
  future_score += calc_future_score_lambda(twobreak.get_vertex(2), twobreak.get_mcolor(), 
                  twobreak.get_vertex(3), twobreak.get_vertex(0));
  future_score += calc_future_score_lambda(twobreak.get_vertex(3), twobreak.get_mcolor(), 
                  twobreak.get_vertex(2), twobreak.get_vertex(1));

  //std::cerr << before_score << " " << future_score << std::endl;
  return std::make_pair(before_score, future_score);         
}


template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_contain_T_consistent_color(vertex_t const & u, vertex_t const & v) const { 
  auto edge_colors = multicolors.split_color_on_tc_color(this->get_all_multicolor_edge(u, v), 20); //graph.count_local_graphs());
  for (auto const & color : edge_colors) { 
    if (multicolors.is_T_consistent_color(color) && !multicolors.is_vec_T_consistent_color(color)) { 
      return true;
    }
  }
  return false; 
} 

template<class mcolor_t>
int BreakpointGraph<mcolor_t>::calculate_cost(vertex_t const & u, vertex_t const & v) const {
  if (u == Infty) {
    mularcs_t mularcs_v = this->get_all_adjacent_vtc_multiedges(v);
    mularcs_v.erase(u);
    return mularcs_v.size();
  }

  if (v == Infty) {
    mularcs_t mularcs_u = this->get_all_adjacent_vtc_multiedges(u);
    mularcs_u.erase(v);
    return mularcs_u.size();
  }

  mularcs_t mularcs_v = this->get_all_adjacent_vtc_multiedges(v);
  mularcs_v.erase(u);
  mularcs_t mularcs_u = this->get_all_adjacent_vtc_multiedges(u);
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

  int count_U = 0; 
  for (auto const & color_set : classes) { 
    count_U += (int) (color_set.second.size() - 1);
  }

  return count_U;
}

template<class mcolor_t>
size_t BreakpointGraph<mcolor_t>::mobility_score_relative_vertex(vertex_t const & x, mcolor_t const & color, vertex_t const & y, vertex_t const & about) const { 
  size_t answer = 0;

  if (this->postponed_deletions.defined(x, y) || !multicolors.is_vec_T_consistent_color(color)) {
    return answer;
  }

  auto calculate_lambda = [&](vertex_t const & my_about, vertex_t const & remove_vertex) -> size_t {
    size_t result = 0;
    mularcs_t mularcs_about; 
#ifdef PSEUDO_EDGE
    if (this->is_pseudo_edge(my_about, color)) {
      mularcs_about = get_clones_adjacent_multiedges_with_info(my_about, color); 
    } else if (my_about != Infty) { 
      mularcs_about = get_real_adjacent_multiedges_with_info(my_about);
    } 
#else
  if (my_about != Infty) { 
    mularcs_about = get_all_adjacent_multiedges_with_info(my_about);
  }
#endif

    mularcs_about.erase(remove_vertex); 
    
    for (auto const & arc : mularcs_about) { 
      if (this->canformQ(arc, color)) {
        ++result;
      }
    } 

    return result;
  }; 

  if (x == about) {  
    answer = calculate_lambda(about, y);
  } else if (y == about) { 
    answer = calculate_lambda(about, x);
  } else { 
    assert((x == about) || (y == about));
  }

  return answer;
}

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const {
  size_t answer = mobility_score_relative_vertex(x, color, y, x) + mobility_score_relative_vertex(x, color, y, y);
  return (answer != 0);    
}

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_mobility_edge(vertex_t const & x, vertex_t const & y) const {
  if ((x == Infty || y == Infty) && is_mobile_irregular_edge) { 
    return true;
  }

  bool non_mobile = false; 
  auto const & central_colors = this->get_all_multicolor_edge_with_info(x, y);
  for (auto color = central_colors.cbegin(); (color != central_colors.cend()) && !non_mobile; ++color) { 
    if (!multicolors.is_T_consistent_color(*color)) { 
      non_mobile = true;
    } else if (!multicolors.is_vec_T_consistent_color(*color)) { 
      non_mobile = true;
    } 
  }

  if (non_mobile) { 
    return false; 
  }

  size_t answer = 0;
  for (auto color = central_colors.cbegin(); color != central_colors.cend(); ++color) { 
    answer += (mobility_score_relative_vertex(x, *color, y, x) + mobility_score_relative_vertex(x, *color, y, y));
  } 

  return (answer != 0);
} 

/*
 *
 *
 */
 #ifdef PSEUDO_EDGE
template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_pseudo_edge(vertex_t const & v, mcolor_t const & color) const { 
  bool answer = false;
  if (mother_verteces.count(v) != 0) {
    auto const & clones = this->mother_verteces.find(v)->second;
    for (auto clone_ind = clones.crbegin(); clone_ind != clones.crend() && !answer; ++clone_ind) { 
      mcolor_t const & clone_color = history.get_clone(*clone_ind).get_mcolor();
      mcolor_t inter_color(clone_color, color, mcolor_t::Intersection);
      if (inter_color.size() > 0 && inter_color.size() < clone_color.size()) { 
        answer = true; 
      }
    } 
  }
  return answer;
}

template<class mcolor_t>
bool BreakpointGraph<mcolor_t>::is_pseudo_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const {
  return (is_pseudo_edge(x, color) || is_pseudo_edge(y, color));
}
#endif 

#endif