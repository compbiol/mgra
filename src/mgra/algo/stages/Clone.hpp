#ifndef CLONE_STAGE_HPP
#define CLONE_STAGE_HPP

namespace algo { 

template<class graph_pack_t>
struct ProcessClone : public algo::AbsStage<graph_pack_t> {  
  using mcolor_t = typename graph_pack_t::mcolor_t;
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using edge_t = typename graph_pack_t::edge_t;  
  using arc_t = typename graph_pack_t::arc_t; 
  using clone_t = typename graph_pack_t::clone_t;
  using ind_arcs_t = typename std::pair<std::pair<vertex_t, mcolor_t>, size_t>;
  
  explicit ProcessClone(size_t max_round)
  : AbsStage<graph_pack_t>("Process clone situation", "clone", round_stage_t, max_round) 
  {
  }
  
  bool run(graph_pack_t & graph_pack) override;

private: 
  using set_arc_t = std::set<arc_t>;

  bool is_good_clones(graph_pack_t const & graph_pack, std::vector<clone_t> const & clones) const;

  std::pair<bool, clone_t> create_clone(graph_pack_t & graph_pack,
    mularcs_t const & mularcs_mother, 
    vertex_t const & v, vertex_t const & father, mularcs_t const & mularcs_father, 
    set_arc_t const & mobile_mother, set_arc_t const & mobile_father);

  std::map<ind_arcs_t, std::set<ind_arcs_t> > split_by_colors(mularcs_t const & mularcs_x, 
      mularcs_t const & mularcs_y) const { 
    utility::equivalence<ind_arcs_t> equiv; 
    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) {
        mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
        if (color.size() > 0) {   
          equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
        } 
      } 
    }  
    equiv.update();
    return equiv.template get_eclasses<std::set<ind_arcs_t> >();    
  }

private:
  DECL_LOGGER("ProcessClone");
};

template<class graph_pack_t>
bool ProcessClone<graph_pack_t>::run(graph_pack_t & graph_pack) {
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 
  
  do {
    number_rear = 0; 
  
    for (vertex_t const & x : graph_pack.graph) {  
      mularcs_t const & mularcs = graph_pack.get_all_adjacent_multiedges(x);

      if (graph_pack.is_duplication_vertex(x) || (mularcs.begin()->second == graph_pack.multicolors.get_complete_color())) {
        continue;
      }  

      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y == Infty || graph_pack.is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!graph_pack.is_mobility_edge(x, y)) { 

          mularcs_t mularcs_x = graph_pack.get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          mularcs_t mularcs_y = graph_pack.get_all_adjacent_multiedges_with_info(y);    
          mularcs_y.erase(x); 

          /*SPLIT ALL EDGES ON MOBILE AND NON MOBILE*/
          set_arc_t mobile_edges_x; 
          set_arc_t non_mobile_edges_x;
          split_by_mobile_property(graph_pack, x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          set_arc_t mobile_edges_y; 
          set_arc_t non_mobile_edges_y;
          split_by_mobile_property(graph_pack, y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          /*CREATE POSSIBLE MOBILE CLONES*/
          std::map<ind_arcs_t, std::set<ind_arcs_t> > classes = split_by_colors(mularcs_x, mularcs_y);          
          std::vector<clone_t> possible_clones;
          
          bool is_all_forks = true; 
          for (auto const & color_set : classes) { 
            mularcs_t mularcs_left; 
            mularcs_t mularcs_right; 
            
            for (auto const & color : color_set.second) { 
              if (color.second == 0) {
                mularcs_left.insert(color.first.first, color.first.second);
              } else {
                mularcs_right.insert(color.first.first, color.first.second);
              } 
            }

            if (mularcs_left.size() == 1 && mularcs_right.size() == 1) { 
              is_all_forks = false; 
            } else if (mularcs_left.size() == 1 && mularcs_left.cbegin()->second == mularcs_right.union_multicolors()) { 
              auto result = create_clone(graph_pack, mularcs_left, x, y, mularcs_right, mobile_edges_x, mobile_edges_y);
              if (result.first) { 
                possible_clones.push_back(result.second);
              } else { 
                is_all_forks = false; 
              }
            } else if (mularcs_right.size() == 1 && mularcs_right.cbegin()->second == mularcs_left.union_multicolors()) { 
              auto result = create_clone(graph_pack, mularcs_right, y, x, mularcs_left, mobile_edges_y, mobile_edges_x);
              if (result.first) { 
                possible_clones.push_back(result.second);
              } else { 
                is_all_forks = false; 
              }
            } else {
              is_all_forks = false; 
            }
          } 

          /*START TO WORK WITH MAIN ALGORITHM*/
          /*CASE 1: Create complete edge*/
          if (is_all_forks) {
            if (is_good_clones(graph_pack, possible_clones)) {
              for (auto const & clone : possible_clones) {
                found = true;
                graph_pack.apply(clone);
                ++number_rear;
              }

              if (found) {
                assert(graph_pack.get_all_multicolor_edge(x, y).empty() || graph_pack.get_all_multicolor_edge(x, y) == graph_pack.multicolors.get_complete_color()); 
              }
            }
          }

          /*CASE 2: Create min T-consistent color*/          
          if (!found) {
            mcolor_t additional_color = graph_pack.multicolors.get_min_addit_color_for_tc(graph_pack.get_all_multicolor_edge(x, y));
            
            if (!additional_color.empty() && additional_color != graph_pack.multicolors.get_complete_color()) {
              std::vector<clone_t> included_clones;

              for (auto clone = possible_clones.cbegin(); (clone != possible_clones.cend()) && !additional_color.empty(); ++clone) { 
                mcolor_t action_color = clone->get_mcolor();
                if (additional_color.includes(action_color) && graph_pack.multicolors.is_vec_T_consistent_color(action_color)) {
                  included_clones.push_back(*clone);
                  additional_color = mcolor_t(additional_color, action_color, mcolor_t::Difference);
                } 
              }

              if (additional_color.empty() && is_good_clones(graph_pack, included_clones)) { 
                for (auto const & clone : included_clones) {
                  graph_pack.apply(clone);
                  found = true;
                  ++number_rear;
                }

                if (found) { 
                  assert(graph_pack.get_all_multicolor_edge(x, y).empty() || graph_pack.multicolors.is_T_consistent_color(graph_pack.get_all_multicolor_edge(x, y)));
                } 
              } 
            } 
          }
        }
      }
    }
    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 

  return isChanged;
}        

template<class graph_pack_t>
std::pair<bool, typename graph_pack_t::clone_t> ProcessClone<graph_pack_t>::create_clone(graph_pack_t & graph_pack, 
    mularcs_t const & mularcs_mother, 
    vertex_t const & v, vertex_t const & father, mularcs_t const & mularcs_father, 
    set_arc_t const & mobile_mother, set_arc_t const & mobile_father) { 

  bool sligshot = (mularcs_mother.size() == 1) && (mularcs_father.size() != 1) 
      && graph_pack.multicolors.is_vec_T_consistent_color(mularcs_mother.cbegin()->second)
      && (mobile_mother.count(*mularcs_mother.cbegin()) != 0);
  size_t count_mobile_father = 0;
  for (auto arc = mularcs_father.cbegin(); arc != (mularcs_father.cend()) && sligshot; ++arc) {
    sligshot = (graph_pack.multicolors.is_vec_T_consistent_color(arc->second));// && (mobile_father.count(*arc) != 0));
    if (mobile_father.count(*arc) != 0) { 
      ++count_mobile_father;
    } 
  }
  assert(mularcs_mother.cbegin()->second == mularcs_father.union_multicolors());   
  
  clone_t clone;
  bool result = false;
  if (sligshot && count_mobile_father > 0) { 
    vertex_t const & mother = mularcs_mother.begin()->first;        
    if (mother == Infty) { 
      vertex_t pseudo_vertex = graph_pack.request_pseudo_infinity_vertex(); 
      clone = clone_t(edge_t(father, v), mularcs_father, arc_t(pseudo_vertex, mularcs_mother.cbegin()->second), true);  
      result = true;
    } else {
      clone = clone_t(edge_t(father, v), mularcs_father, *(mularcs_mother.cbegin()), false);  
      result = true;
    }  
  }

  return std::make_pair(result, clone);
}


template<class graph_pack_t>
bool ProcessClone<graph_pack_t>::is_good_clones(graph_pack_t const & graph_pack, 
  std::vector<clone_t> const & clones) const {
  bool is_all_good = true;

  for (auto vec_color = graph_pack.multicolors.cbegin_vec_T_consistent_color(); vec_color != graph_pack.multicolors.cend_vec_T_consistent_color(); ++vec_color) {
    size_t count_diff = 0;
    
    mcolor_t vec_target_color = *vec_color;
    for (clone_t const & clone : clones) { 
      mcolor_t local_color = clone.get_mcolor();
      if (vec_target_color.includes(local_color)) {
        vec_target_color = mcolor_t(vec_target_color, local_color, mcolor_t::Difference);
        ++count_diff;
      }
    }

    if (vec_target_color.empty() && (count_diff != 1)) {
      is_all_good = false; 
    }
  }

  return is_all_good;
}

} 

#endif
