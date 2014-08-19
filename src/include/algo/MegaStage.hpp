#ifndef MEGASTAGE_HPP
#define MEGASTAGE_HPP

template<class graph_t>
struct Algorithm<graph_t>::MegaStage : public Algorithm<graph_t>::Stage {
  typedef Stage base;
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename std::pair<std::pair<vertex_t, mcolor_t>, size_t> ind_acrs_t;
            
  typedef typename graph_t::clone_t clone_t;        
  typedef typename graph_t::twobreak_t twobreak_t;


  explicit MegaStage(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  , pseudo_infinity_vertex(0)
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override { 
    return "Mega stage.";
  }

private: 
  typedef std::set<edge_t> set_edge_t;

  void split_by_mobile_property(vertex_t const & v, mularcs_t const & mularcs, 
      set_edge_t& mobiles, set_edge_t& non_mobiles) const {

    for (auto const & arc : mularcs) { 
      if (this->graph->is_mobility_edge(v, arc.second, arc.first)) { 
        mobiles.insert(arc); 
      } else { 
        non_mobiles.insert(arc);
      }
    }   
  }

  std::map<ind_acrs_t, std::set<ind_acrs_t> > split_by_colors(mularcs_t const & mularcs_x, 
      mularcs_t const & mularcs_y) const { 

    utility::equivalence<ind_acrs_t> equiv; 
    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) {
        //if (this->graph->is_vec_T_consistent_color(arc_x.second) && this->graph->is_vec_T_consistent_color(arc_y.second))
        mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
        if (color.size() > 0) {   
          equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
        } 
      } 
    }  
    equiv.update();
    return equiv.template get_eclasses<std::set<ind_acrs_t> >();    
  }

  std::pair<bool, clone_t> create_clone(mularcs_t const & mularcs_mother, vertex_t const & v, 
                  vertex_t const & father, mularcs_t const & mularcs_father,
                  set_edge_t const & mobile_mother, set_edge_t const & mobile_father);

  template<class action_t>
  std::set<mcolor_t> get_worked_colors(std::vector<action_t> actions) { 
    std::set<mcolor_t> result;
    for (auto const & action : actions) { 
      result.insert(action.get_mcolor());
    }
    return result; 
  }

  bool is_good_actions(std::set<mcolor_t> const & actions) const;

  mcolor_t get_min_addit_color_for_tc(mcolor_t const & color) const { 
    mcolor_t min_color = this->graph->get_complete_color();
    for (auto col = this->graph->cbegin_T_consistent_color(); col != this->graph->cend_T_consistent_color(); ++col) {
      if (col->includes(color)) {
        mcolor_t diff_color(*col, color, mcolor_t::Difference);
        if (diff_color.size() < min_color.size()) {
          min_color = diff_color;
        } 
      } 
    } 
    return min_color;
  }
private: 
  size_t pseudo_infinity_vertex;  
};

/*
Idea: 
  If all edges mobile, is it true that all split on clones and verteces
*/
template<class graph_t>
bool Algorithm<graph_t>::MegaStage::do_action() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 

    for (vertex_t const & x: *this->graph) {  
      mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x);

      if (this->graph->is_duplication_vertex(x)) {
        continue;
      }  

      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y == Infty || this->graph->is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!this->graph->is_mobility_edge(x, y) && im->second != this->graph->get_complete_color()) { 
          mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y);    
          mularcs_y.erase(x);  

          /*SPLIT ALL VERTEX ON MOBILE AND NON MOBILE*/
          set_edge_t mobile_edges_x; 
          set_edge_t non_mobile_edges_x;
          split_by_mobile_property(x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          set_edge_t mobile_edges_y; 
          set_edge_t non_mobile_edges_y;
          split_by_mobile_property(y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          /*SPLIT ALL COLORS ON CLONE AND TWO-BREAKS*/
          std::map<ind_acrs_t, std::set<ind_acrs_t> > classes = split_by_colors(mularcs_x, mularcs_y);          

          /*CREATE POSSIBLE MOBILE TWO-BREAKS and CLONES*/
          std::vector<twobreak_t> possible_twobreaks;
          std::vector<clone_t> possible_clones;

          size_t bad_edges = 0;
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

            /*std::cerr << "Start new class" << std::endl;
            std::cerr << mularcs_left.size() << " " << mularcs_right.size() << std::endl;

            std::cerr << "mulacrs_left vertex have " << std::endl; 
            for (auto const & l : mularcs_left) { 
              std::cerr << genome_match::mcolor_to_name(l.second) << " " << this->graph->is_vec_T_consistent_color(l.second) << " " << l.first << std::endl;
            }

            std::cerr << "mularcs_right vertex have " << std::endl; 
            for (auto const & r : mularcs_right) {
              std::cerr << genome_match::mcolor_to_name(r.second) << " " << this->graph->is_vec_T_consistent_color(r.second) << " " << r.first << std::endl;
            }*/

            if (mularcs_left.size() == 1 && mularcs_right.size() == 1 && mularcs_left.cbegin()->second == mularcs_right.cbegin()->second) { 
              mcolor_t const & color = mularcs_left.cbegin()->second;
              if (mobile_edges_x.count(*mularcs_left.cbegin()) != 0 
                && mobile_edges_y.count(*mularcs_right.cbegin()) != 0
                && this->graph->is_vec_T_consistent_color(color)) { 
                vertex_t const & u = mularcs_left.cbegin()->first; 
                vertex_t const & v = mularcs_right.cbegin()->first;
                twobreak_t twobreak(x, u, y, v, color);
                possible_twobreaks.push_back(twobreak);
              } else { 
                ++bad_edges;
              }
            } else if (mularcs_left.size() == 1) { 
              auto result = create_clone(mularcs_left, x, y, mularcs_right, mobile_edges_x, mobile_edges_y);
              if (result.first) { 
                possible_clones.push_back(result.second);
              } else { 
                ++bad_edges;
              }
            } else if (mularcs_right.size() == 1) { 
              auto result = create_clone(mularcs_right, y, x, mularcs_left, mobile_edges_y, mobile_edges_x);
              if (result.first) { 
                possible_clones.push_back(result.second);
              } else { 
                ++bad_edges;
              }
            } else {
              ++bad_edges; 
              //assert(false);
            }
          } 

          /*START TO WORK WITH MAIN ALGORITHM*/
          /*CASE 1: Create complete edge*/
          if (non_mobile_edges_x.empty() && non_mobile_edges_y.empty()) {
            assert(bad_edges == 0);
            std::set<mcolor_t> actions = get_worked_colors(possible_twobreaks);
            std::set<mcolor_t> clone_action = get_worked_colors(possible_clones);
            actions.insert(clone_action.begin(), clone_action.end()); 

            if (is_good_actions(actions)) { 
              //std::cerr << "Do cloning and two-break in first case" << std::endl;
              for (auto const & twobreak : possible_twobreaks) {
                //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                this->graph->apply(twobreak);
                found = true;
                ++number_rear;
              } 
              
              for (auto const & clone : possible_clones) {
                this->graph->apply(clone);
                found = true;
                ++number_rear;
              }
              //std::cerr << genome_match::mcolor_to_name(this->graph->get_all_multicolor_edge(x, y)) << std::endl;
              assert(this->graph->get_all_multicolor_edge(x, y) == this->graph->get_complete_color());
            }
          } 

          /*CASE 2: Create min T-consistent color*/          
          if (!found) { 
            mcolor_t additional_color = get_min_addit_color_for_tc(im->second);
            if (!additional_color.empty() && additional_color != this->graph->get_complete_color()) {
              std::vector<twobreak_t> included_twobreaks;
              std::vector<clone_t> included_clones;
              
              for (auto twobreak = possible_twobreaks.cbegin(); (twobreak != possible_twobreaks.cend()) && !additional_color.empty(); ++twobreak) { 
                mcolor_t action_color = twobreak->get_mcolor();
                if (additional_color.includes(action_color) && this->graph->is_vec_T_consistent_color(action_color)) {
                  included_twobreaks.push_back(*twobreak);
                  additional_color = mcolor_t(additional_color, action_color, mcolor_t::Difference);
                } 
              }

              /*for (auto clone = possible_clones.cbegin(); (clone != possible_clones.cend()) && !additional_color.empty(); ++clone) { 
                mcolor_t action_color = clone->get_mcolor();
                if (additional_color.includes(action_color) && this->graph->is_vec_T_consistent_color(action_color)) {
                  included_clones.push_back(*clone);
                  additional_color = mcolor_t(additional_color, action_color, mcolor_t::Difference);
                }
              }*/

              if (additional_color.empty()) {
                std::set<mcolor_t> actions = get_worked_colors(included_twobreaks);
                std::set<mcolor_t> clone_action = get_worked_colors(included_clones);
                actions.insert(clone_action.begin(), clone_action.end()); 
             
                if (is_good_actions(actions)) {                
                  //std::cerr << "Do cloning and two-break in second case" << std::endl;
                  for (auto const & twobreak : included_twobreaks) {
                    //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                    this->graph->apply(twobreak);
                    found = true;
                    ++number_rear;
                  } 
              
                  /*for (auto const & clone : included_clones) {
                    this->graph->apply(clone);
                    found = true;
                    ++number_rear;
                  }*/
                  
                  //assert(this->graph->is_T_consistent_color(this->graph->get_all_multicolor_edge(x, y)));
                }
              }
            } 
          }

          /*CASE 3: Process unique two-breaks*/          
          if (!found && !this->graph->is_vec_T_consistent_color(im->second)) {
            /*auto const count_variant_Lambda = [&] (arc_t const & viewed, mcolor_t const & Q, arc_t const & remove) -> size_t { 
              size_t number_variant = 0; 
              if (viewed.first != Infty) {
                mularcs_t&& end_f = this->graph->get_all_adjacent_multiedges_with_info(viewed.first);
                end_f.erase(viewed.second);
                end_f.erase(remove.first); 
                end_f.erase(remove.second); 
                for (auto arc = end_f.cbegin(); arc != end_f.cend(); ++arc) { 
                  if (this->graph->canformQ(*arc, Q)) {
                    ++number_variant;
                  } 
                } 
              } 

              if (viewed.second != Infty) {                
                mularcs_t&& end_s = this->graph->get_all_adjacent_multiedges_with_info(viewed.second);
                end_s.erase(viewed.first);
                end_s.erase(remove.first); 
                end_s.erase(remove.second);
                for (auto arc = end_s.cbegin(); arc != end_s.cend(); ++arc) { 
                  if (this->graph->canformQ(*arc, Q)) {
                    ++number_variant;
                  }
                }   
              } 

              return number_variant; 
            };*/

            for (auto twobreak = possible_twobreaks.cbegin(); (twobreak != possible_twobreaks.cend()) && !found; ++twobreak) { 
              size_t count = 42; //count_variant_Lambda(twobreak->get_arc(0), twobreak->get_mcolor(), twobreak->get_arc(1));
              //count += count_variant_Lambda(twobreak->get_arc(1), twobreak->get_mcolor(), twobreak->get_arc(0));    
              if (count == 0) {
                  bool flag = true; 

                  auto temp_color = twobreak->get_mcolor();
                  if (this->graph->is_pseudo_edge(twobreak->get_vertex(0), temp_color) || this->graph->is_pseudo_edge(twobreak->get_vertex(3), temp_color)) { 
                    flag = (this->graph->is_pseudo_edge(twobreak->get_vertex(0), temp_color) && 
                          this->graph->get_clones_adjacent_multiedges_with_info(twobreak->get_vertex(0), temp_color).defined(twobreak->get_vertex(3))) 
                      || (this->graph->is_pseudo_edge(twobreak->get_vertex(3), temp_color) 
                        && this->graph->get_clones_adjacent_multiedges_with_info(twobreak->get_vertex(3), temp_color).defined(twobreak->get_vertex(0)));
                  } 

                if (flag) { 

                //std::cerr << " Do 2-break in third case " << std::endl;
                //std::cerr << twobreak->get_vertex(0) << " " << twobreak->get_vertex(1) << " " << twobreak->get_vertex(2) << " " << twobreak->get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak->get_mcolor()) << std::endl;
                this->graph->apply(*twobreak);  
                found = true;
                ++number_rear;
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

  return (number_rear != 0);
}

template<class graph_t>
std::pair<bool, typename graph_t::clone_t> Algorithm<graph_t>::MegaStage::create_clone(mularcs_t const & mularcs_mother, 
    vertex_t const & v, vertex_t const & father, mularcs_t const & mularcs_father, 
    set_edge_t const & mobile_mother, set_edge_t const & mobile_father) { 

  bool sligshot = (mularcs_mother.size() == 1) && (mularcs_father.size() != 1) 
      && this->graph->is_vec_T_consistent_color(mularcs_mother.cbegin()->second)
      && (mobile_mother.count(*mularcs_mother.cbegin()) != 0);
  for (auto arc = mularcs_father.cbegin(); arc != (mularcs_father.cend()) && sligshot; ++arc) {
    sligshot = (this->graph->is_vec_T_consistent_color(arc->second) && (mobile_father.count(*arc) != 0));
  }
  assert(mularcs_mother.cbegin()->second == mularcs_father.union_multicolors());   
  
  clone_t clone;
  bool result = false;
  if (sligshot) { 
    vertex_t const & mother = mularcs_mother.begin()->first;        
    if (mother == Infty) { 
      std::string pseudo_vertex = "o" + std::to_string(pseudo_infinity_vertex) + "o";
      ++pseudo_infinity_vertex;
      clone = clone_t(arc_t(father, v), mularcs_father, edge_t(pseudo_vertex, mularcs_mother.cbegin()->second), true);  
      result = true;
    } else {
      clone = clone_t(arc_t(father, v), mularcs_father, *(mularcs_mother.cbegin()), false);  
      result = true;
    }  
  }

  return std::make_pair(result, clone);
}

template<class graph_t>
bool Algorithm<graph_t>::MegaStage::is_good_actions(std::set<mcolor_t> const & actions) const {
  bool is_all_good = true;

  for (auto vec_color = this->graph->cbegin_vec_T_consistent_color(); vec_color != this->graph->cend_vec_T_consistent_color(); ++vec_color) {
    size_t count_diff = 0;
    mcolor_t vec_target_color = *vec_color;
    for (mcolor_t const & local_color : actions) { 
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

#endif