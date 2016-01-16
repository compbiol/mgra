#ifndef PROPERTY_FUNC_IMPL_HPP
#define PROPERTY_FUNC_IMPL_HPP

/*
 * Implementation of the function about property for vertex
 */
template<class mcolor_t>
bool GraphPack<mcolor_t>::is_simple_vertex(vertex_t const & v) const {
  return (this->graph.degree_vertex(v) == 2); //(!is_duplication_vertex(v) && !is_indel_vertex(v) &&
}  

template<class mcolor_t>
bool GraphPack<mcolor_t>::is_have_self_loop(vertex_t const & v) const {
  return !(this->get_all_multicolor_edge(v, v).empty()); 
} 
 
template<class mcolor_t>
bool GraphPack<mcolor_t>::is_indel_vertex(vertex_t const & v) const {
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
bool GraphPack<mcolor_t>::is_duplication_vertex(vertex_t const & v) const {  
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
 * Implementation of the function for calculate different scores. 
 * Mobility score, verteces score, twobreak score. 
 * Detailed about this scores see in journal paper corresponding MGRA 2.0.
 */
template<class mcolor_t>
size_t GraphPack<mcolor_t>::mobility_score(edge_t const & viewed, mcolor_t const & color, edge_t const & removed) const {
  auto calculate_lambda = [&] (vertex_t const & target, vertex_t const & other) -> size_t {
    mularcs_t mularcs; 
    
    if (target != Infty) { 
      mularcs = this->get_all_adjacent_multiedges_with_info(target); 
    } 

    mularcs.erase(other);
    mularcs.erase(removed.first); 
    mularcs.erase(removed.second); 
    
    size_t number_variant = 0; 
    for (auto const & arc : mularcs) { 
      if (this->canformQ(arc.first, color)) {
        ++number_variant;
      } 
    }
    return number_variant;
  };

  return (calculate_lambda(viewed.first, viewed.second) + calculate_lambda(viewed.second, viewed.first)); 
}

template<class mcolor_t>
std::pair<size_t, size_t> GraphPack<mcolor_t>::is_decrease_verteces_score(twobreak_t const & twobreak) const { 
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
      mularcs.erase(removed, color); 

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

      return calc_score_lambda(mularcs);
    } else { 
      return 0; 
    }
  };

  size_t before_score = calc_vertex_score_lambda(twobreak.get_vertex(0)) + calc_vertex_score_lambda(twobreak.get_vertex(1))
        + calc_vertex_score_lambda(twobreak.get_vertex(2)) + calc_vertex_score_lambda(twobreak.get_vertex(3));

  size_t future_score = calc_future_score_lambda(twobreak.get_vertex(0), twobreak.get_mcolor(), 
                  twobreak.get_vertex(1), twobreak.get_vertex(2));
  future_score += calc_future_score_lambda(twobreak.get_vertex(1), twobreak.get_mcolor(), 
                  twobreak.get_vertex(0), twobreak.get_vertex(3));
  future_score += calc_future_score_lambda(twobreak.get_vertex(2), twobreak.get_mcolor(), 
                  twobreak.get_vertex(3), twobreak.get_vertex(0));
  future_score += calc_future_score_lambda(twobreak.get_vertex(3), twobreak.get_mcolor(), 
                  twobreak.get_vertex(2), twobreak.get_vertex(1));

  return std::make_pair(before_score, future_score);         
}

template<class mcolor_t>
int GraphPack<mcolor_t>::calculate_cost(vertex_t const & u, vertex_t const & v) const {
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
    
  using ind_acr_t = std::tuple<vertex_t, mcolor_t, size_t>; 
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
size_t GraphPack<mcolor_t>::mobility_score_relative_vertex(vertex_t const & x, mcolor_t const & color, vertex_t const & y, vertex_t const & about) const { 
  size_t answer = 0;

  if (!multicolors.is_vec_T_consistent_color(color)) {
    return answer;
  }

  auto calculate_lambda = [&](vertex_t const & my_about, vertex_t const & remove_vertex) -> size_t {
    size_t result = 0;
    mularcs_t mularcs_about; 
    
    if (my_about != Infty) { 
      mularcs_about = this->get_all_adjacent_multiedges_with_info(my_about);
    }

    mularcs_about.erase(remove_vertex); 
    
    for (auto const & arc : mularcs_about) { 
      if (this->canformQ(arc.first, color)) {
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

/*
 * Implementation of the function for different edges property. 
 */
template<class mcolor_t>
bool GraphPack<mcolor_t>::is_contain_T_consistent_color(vertex_t const & u, vertex_t const & v) const { 
  auto edge_colors = multicolors.split_color_on_tc_color(this->get_all_multicolor_edge(u, v), graph.count_local_graphs() + 1); 
  for (auto const & color : edge_colors) { 
    if (multicolors.is_T_consistent_color(color) && !multicolors.is_vec_T_consistent_color(color)) { 
      return true;
    }
  }
  return false; 
} 

template<class mcolor_t>
bool GraphPack<mcolor_t>::is_mobility_edge(vertex_t const & x, mcolor_t const & color, vertex_t const & y) const {
  size_t answer = mobility_score_relative_vertex(x, color, y, x) + mobility_score_relative_vertex(x, color, y, y);
  return (answer != 0);    
}

template<class mcolor_t>
bool GraphPack<mcolor_t>::is_mobility_edge(vertex_t const & x, vertex_t const & y) const {
  if (x == Infty || y == Infty) { 
    return !is_contain_T_consistent_color(x, y);
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

/**
 * Fixme need to think about it.
 */
template<class mcolor_t>
bool GraphPack<mcolor_t>::is_consistency_graph() const {
  return (mother_verteces.size() == 0);
}

/*
 * Implementation of the function for say canform or not. 
 * "canform" definition see in journal paper corresponding MGRA 2.0.
 */

/*
 * Method is return a series of two-breaks, which needs for create edge (x, y)
 * with multicolor = COLOR 
 */
template<class mcolor_t>
typename GraphPack<mcolor_t>::transform_t GraphPack<mcolor_t>::take_edge_on_color(
              vertex_t const & x, mcolor_t const & color, vertex_t const & y) const { 
  std::list<twobreak_t> result;

  if (y == Infty) {
    mcolor_t need_color(color, get_all_multicolor_edge(x, y), mcolor_t::Difference);    
    mularcs_t mularcs_x = get_all_adjacent_multiedges_with_info(x);
    mularcs_x.erase(y);
  
    for (arc_t const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && multicolors.is_vec_T_consistent_color(arc.second)) {
        twobreak_t twobreak(x, arc.first, Infty, Infty, arc.second);
        result.push_back(twobreak);
      } 
    }
  } else { 
    mcolor_t need_color(color, get_all_multicolor_edge(x, y), mcolor_t::Difference);
    
    mularcs_t mularcs_x = get_all_adjacent_multiedges_with_info(x);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = get_all_adjacent_multiedges_with_info(y);
    mularcs_y.erase(x);

    using colacr_t = std::pair<std::pair<vertex_t, structure::Mcolor>, size_t>;

    utility::equivalence<colacr_t> equiv; 
    for (arc_t const & arc_x : mularcs_x) { 
      for (arc_t const & arc_y : mularcs_y) {
        if (need_color.includes(arc_x.second) && need_color.includes(arc_y.second)) {  
          mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
          if (color.size() > 0) {   
            equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
          }
        }
      } 
    }      
    equiv.update();
    
    std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.get_eclasses<std::set<colacr_t> >();   
   
    mcolor_t temp_need_color = need_color;
    for (auto const & color_set : classes) {
      mularcs_t left; 
      mularcs_t right; 

      for (auto const & color : color_set.second) { 
        if (color.second == 0) {
          left.insert(color.first.first, color.first.second);
        } else {
          right.insert(color.first.first, color.first.second);
        }
      }

      if (left.size() == 1) {
        temp_need_color = mcolor_t(temp_need_color, left.begin()->second, mcolor_t::Difference); 
      } else if (right.size() == 1) { 
        temp_need_color = mcolor_t(temp_need_color, right.begin()->second, mcolor_t::Difference); 
      } 
    } 

    if (!temp_need_color.empty()) { 
      return result; 
    }

    temp_need_color = need_color;
    for (auto const & color_set : classes) { 
      std::multimap<mcolor_t, vertex_t> left; 
      std::multimap<mcolor_t, vertex_t> right; 

      for (auto const & color : color_set.second) { 
        if (color.second == 0) {
          left.insert(std::make_pair(color.first.second, color.first.first));
        } else {
          right.insert(std::make_pair(color.first.second, color.first.first));
        } 
      }
    
      if (left.size() == 1 && right.size() == 1 && left.begin()->first == right.begin()->first) {
        assert(multicolors.is_vec_T_consistent_color(left.begin()->first)); 
        assert(multicolors.is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          twobreak_t twobreak(x, left.begin()->second, y, right.begin()->second, left.begin()->first);
          temp_need_color = mcolor_t(temp_need_color, left.begin()->first, mcolor_t::Difference);
          result.push_back(twobreak);
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          assert(multicolors.is_vec_T_consistent_color(left.begin()->first)); 
          
          std::list<twobreak_t> local_twobreaks; 

          auto iter_right = right.begin();
          while ((local_twobreaks.empty()) && (iter_right != right.end())) { 
            local_twobreaks = take_edge_on_color(y, left.begin()->first, iter_right->second);
            if (local_twobreaks.empty()) { 
              ++iter_right;
            } 
          } 
          
          if (iter_right != right.end()) {
            result.insert(result.end(), local_twobreaks.begin(), local_twobreaks.end());
            twobreak_t twobreak(x, left.begin()->second, y, iter_right->second, left.begin()->first);
            temp_need_color = mcolor_t(temp_need_color, left.begin()->first, mcolor_t::Difference);
            result.push_back(twobreak);
          }          
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          assert(multicolors.is_vec_T_consistent_color(right.begin()->first));  
          
          std::list<twobreak_t> local_twobreaks;
          auto iter_left = left.begin();
          while ((local_twobreaks.empty()) && (iter_left != left.end())) {
            local_twobreaks = take_edge_on_color(x, right.begin()->first, iter_left->second);
            if (local_twobreaks.empty()) { 
              ++iter_left;
            } 
          }

          if (iter_left != left.end()) {
            result.insert(result.end(), local_twobreaks.begin(), local_twobreaks.end());
            twobreak_t twobreak(x, iter_left->second, y, right.begin()->second, right.begin()->first);
            temp_need_color = mcolor_t(temp_need_color, right.begin()->first, mcolor_t::Difference);
            result.push_back(twobreak);
          }             
        } else { 
          ;
        }   
      } else { 
        assert(left.size() == 1 || right.size() == 1); 
      } 
    }

    if (!temp_need_color.empty()) { 
      result.clear(); 
    } 
  } 

  return result; 
}

#endif