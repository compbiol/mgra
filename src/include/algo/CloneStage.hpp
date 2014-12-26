#ifndef CLONE_STAGE_HPP
#define CLONE_STAGE_HPP

template<class graph_t>
struct Algorithm<graph_t>::ProcessClone : public Algorithm<graph_t>::Stage {
  typedef Stage base;
  
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::edge_t edge_t;  
  typedef typename graph_t::arc_t arc_t;  
  typedef typename graph_t::clone_t clone_t;

  typedef typename std::pair<std::pair<vertex_t, mcolor_t>, size_t> ind_acrs_t;
  
  explicit ProcessClone(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override { 
    return "Process clone situation.";
  }

private: 
  typedef std::set<arc_t> set_arc_t;

  bool is_good_clones(std::vector<clone_t> const & clones) const;

  std::pair<bool, clone_t> create_clone(mularcs_t const & mularcs_mother, 
    vertex_t const & v, vertex_t const & father, mularcs_t const & mularcs_father, 
    set_arc_t const & mobile_mother, set_arc_t const & mobile_father);

  std::map<ind_acrs_t, std::set<ind_acrs_t> > split_by_colors(mularcs_t const & mularcs_x, 
      mularcs_t const & mularcs_y) const { 
    utility::equivalence<ind_acrs_t> equiv; 
    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) {
        mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
        if (color.size() > 0) {   
          equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
        } 
      } 
    }  
    equiv.update();
    return equiv.template get_eclasses<std::set<ind_acrs_t> >();    
  }

private:
  DECL_LOGGER("ProcessClone");
};

template<class graph_t>
bool Algorithm<graph_t>::ProcessClone::do_action() {
  INFO("Process clone situation in fair edge.") 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 
  
  do {
    number_rear = 0; 
  
    for(vertex_t const & x : *this->graph) {  
      mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x);

      if (this->graph->is_duplication_vertex(x) || (mularcs.begin()->second == this->graph->get_complete_color())) {
        continue;
      }  

      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y == Infty || this->graph->is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!this->graph->is_mobility_edge(x, y)) { 

          mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y);    
          mularcs_y.erase(x); 

          /*SPLIT ALL EDGES ON MOBILE AND NON MOBILE*/
          set_arc_t mobile_edges_x; 
          set_arc_t non_mobile_edges_x;
          this->split_by_mobile_property(x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          set_arc_t mobile_edges_y; 
          set_arc_t non_mobile_edges_y;
          this->split_by_mobile_property(y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          /*CREATE POSSIBLE MOBILE CLONES*/
          std::map<ind_acrs_t, std::set<ind_acrs_t> > classes = split_by_colors(mularcs_x, mularcs_y);          
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
              auto result = create_clone(mularcs_left, x, y, mularcs_right, mobile_edges_x, mobile_edges_y);
              if (result.first) { 
                possible_clones.push_back(result.second);
              } else { 
                is_all_forks = false; 
              }
            } else if (mularcs_right.size() == 1 && mularcs_right.cbegin()->second == mularcs_left.union_multicolors()) { 
              auto result = create_clone(mularcs_right, y, x, mularcs_left, mobile_edges_y, mobile_edges_x);
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
            if (is_good_clones(possible_clones)) {
              for (auto const & clone : possible_clones) {
                //std::cerr << "Do clone in first case" << std::endl;
                found = true;
                this->graph->apply(clone);
                ++number_rear;
              }

              if (found) {
                assert(this->graph->get_all_multicolor_edge(x, y).empty() || this->graph->get_all_multicolor_edge(x, y) == this->graph->get_complete_color()); 
              }
            }
          }

          /*CASE 2: Create min T-consistent color*/          
          if (!found) {
            mcolor_t additional_color = this->get_min_addit_color_for_tc(this->graph->get_all_multicolor_edge(x, y));
            
            if (!additional_color.empty() && additional_color != this->graph->get_complete_color()) {
              std::vector<clone_t> included_clones;

              for (auto clone = possible_clones.cbegin(); (clone != possible_clones.cend()) && !additional_color.empty(); ++clone) { 
                mcolor_t action_color = clone->get_mcolor();
                if (additional_color.includes(action_color) && this->graph->is_vec_T_consistent_color(action_color)) {
                  included_clones.push_back(*clone);
                  additional_color = mcolor_t(additional_color, action_color, mcolor_t::Difference);
                } 
              }

              if (additional_color.empty() && is_good_clones(included_clones)) { 
                for (auto const & clone : included_clones) {
                  //std::cerr << " Do 2-break in second case " << std::endl;
                  //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                  this->graph->apply(clone);
                  found = true;
                  ++number_rear;
                }

                if (found) { 
                  assert(this->graph->get_all_multicolor_edge(x, y).empty() || this->graph->is_T_consistent_color(this->graph->get_all_multicolor_edge(x, y)));
                } 
              } 
            } 
          } 

          /*CASE 3: Mobility of one*/
          /*if (!found) {
            for (auto const & clone : possible_clones) {
              size_t count = this->graph->mobility_score(edge_t(clone.get_central_arc().second, clone.get_mother_edge().first), clone.get_mcolor(), edge_t(clone.get_central_arc().first, "")); 
              
              if (count == 0) { 
                found = true;
                this->graph->apply(clone);
                ++number_rear;
              } 
            } 
          }*/

        }
      }
    }
    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 

  return isChanged;
}        

template<class graph_t>
std::pair<bool, typename graph_t::clone_t> Algorithm<graph_t>::ProcessClone::create_clone(mularcs_t const & mularcs_mother, 
    vertex_t const & v, vertex_t const & father, mularcs_t const & mularcs_father, 
    set_arc_t const & mobile_mother, set_arc_t const & mobile_father) { 

  bool sligshot = (mularcs_mother.size() == 1) && (mularcs_father.size() != 1) 
      && this->graph->is_vec_T_consistent_color(mularcs_mother.cbegin()->second)
      && (mobile_mother.count(*mularcs_mother.cbegin()) != 0);
  size_t count_mobile_father = 0;
  for (auto arc = mularcs_father.cbegin(); arc != (mularcs_father.cend()) && sligshot; ++arc) {
    sligshot = (this->graph->is_vec_T_consistent_color(arc->second));// && (mobile_father.count(*arc) != 0));
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
      vertex_t pseudo_vertex = "o" + std::to_string(this->pseudo_infinity_vertex) + "o";
      ++this->pseudo_infinity_vertex;
      clone = clone_t(edge_t(father, v), mularcs_father, arc_t(pseudo_vertex, mularcs_mother.cbegin()->second), true);  
      result = true;
    } else {
      clone = clone_t(edge_t(father, v), mularcs_father, *(mularcs_mother.cbegin()), false);  
      result = true;
    }  
  }

  return std::make_pair(result, clone);
}


template<class graph_t>
bool Algorithm<graph_t>::ProcessClone::is_good_clones(std::vector<clone_t> const & clones) const {
  bool is_all_good = true;

  for (auto vec_color = this->graph->cbegin_vec_T_consistent_color(); vec_color != this->graph->cend_vec_T_consistent_color(); ++vec_color) {
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
      //(!this->graph->canformQ(clones.begin()->get_central_arc().first, *vec_color) 
      //  || !this->graph->canformQ(clones.begin()->get_central_arc().second, *vec_color));
    }
  }

  return is_all_good;
}

#endif
