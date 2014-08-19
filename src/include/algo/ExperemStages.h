#ifndef EXPEREM_STAGE_H_
#define EXPEREM_STAGE_H_

template<class graph_t>
struct Algorithm<graph_t>::ExperemProcessClone : public Algorithm<graph_t>::Stage {
  typedef Stage base;
  
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::clone_t clone_t;        
  
  typedef typename std::pair<std::pair<vertex_t, mcolor_t>, size_t> ind_acrs_t;
  

  explicit ExperemProcessClone(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  , pseudo_infinity_vertex(0)
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override { 
    return "Experem process clone situation.";
  }

private: 
  typedef std::set<edge_t> set_edge_t;

  bool is_good_clones(std::vector<clone_t> const & clones) const;

  std::pair<bool, clone_t> create_clone(mularcs_t const & mularcs_mother, 
    vertex_t const & v, vertex_t const & father, mularcs_t const & mularcs_father, 
    set_edge_t const & mobile_mother, set_edge_t const & mobile_father);

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

private: 
  size_t pseudo_infinity_vertex;  
};

template<class graph_t>
bool Algorithm<graph_t>::ExperemProcessClone::do_action() { 
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
          set_edge_t mobile_edges_x; 
          set_edge_t non_mobile_edges_x;
          split_by_mobile_property(x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          set_edge_t mobile_edges_y; 
          set_edge_t non_mobile_edges_y;
          split_by_mobile_property(y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          /*CREATE POSSIBLE MOBILE CLONES*/
          /*SPLIT ALL COLORS ON CLONE AND TWO-BREAKS*/
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
          if (non_mobile_edges_x.empty() && non_mobile_edges_y.empty() && is_all_forks) {
            if (is_good_clones(possible_clones)) {
              for (auto const & clone : possible_clones) {
                //std::cerr << "Do clone in first case" << std::endl;
                found = true;
                this->graph->apply(clone);
                ++number_rear;
              } 
            }
          }

          /*CASE 2: Create edge with tip*/          
          /*if (!found && !non_mobile_edges_x.empty() && !non_mobile_edges_y.empty()) {
            ;
            for (auto edge = non_mobile_edges_x.begin(); edge != non_mobile_edges_x.end() && !found; ++edge) { 
              if (this->graph->degree_vertex(edge->first) == 2) {
                mularcs_t mularcs_p = this->graph->get_all_adjacent_multiedges(edge->first); 
                mularcs_p.erase(x);

                if (this->graph->is_vec_T_consistent_color(mularcs_p.begin()->second) && 
                  non_mobile_edges_y.count(edge_t(mularcs_p.begin()->first, edge->second)) != 0) { 
                  mcolor_t target_color = mularcs_p.begin()->second;

                  std::cerr << "central edge " << x << " " << y << std::endl;
                  std::cerr << "Edge (p, q): " << edge->first << " " << mularcs_p.begin()->first << " " << genome_match::mcolor_to_name(target_color) << std::endl;
                  std::cerr << "Edge (x, p) and (y, q) have color: " << genome_match::mcolor_to_name(edge->second) << std::endl;
                  
                  for (auto const & clone : possible_clones) { 
                    mcolor_t action_color = clone.get_mcolor(); 
                    if (target_color.includes(action_color)) { 
                      std::cerr << " Do clone in fourth case " << std::endl;
                      std::cerr << " Clone color " << genome_match::mcolor_to_name(action_color) << std::endl;
                      this->graph->apply(clone);
                      ++number_rear;
                      found = true; 
                    } 
                  }
                }
              }
            } 
          }*/

          /*CASE 2: Mobility of one*/
          /*if (!found) {
            for (auto const & clone : possible_clones) {
              size_t count = this->graph->mobility_score(arc_t(clone.get_central_arc().second, clone.get_mother_edge().first), clone.get_mcolor(), arc_t(clone.get_central_arc().first, "")); 
              
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
std::pair<bool, typename graph_t::clone_t> Algorithm<graph_t>::ExperemProcessClone::create_clone(mularcs_t const & mularcs_mother, 
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
bool Algorithm<graph_t>::ExperemProcessClone::is_good_clones(std::vector<clone_t> const & clones) const {
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
    }
  }

  return is_all_good;
}

/*
template<class graph_t>
bool Algorithm<graph_t>::ExperemProcessClone::do_action() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 

    for (vertex_t const & x: *this->graph) {  
      mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x);
        
      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y == Infty) {
          continue;
        }

        if (!this->graph->is_mobility_edge(x, y) && im->second != this->graph->get_complete_color()) {  
          auto const filter_lambda = [&] (vertex_t const & s, mularcs_t const & mul, std::set<edge_t>& mobiles, std::set<edge_t>& non_mobiles) -> void {
            for (auto arc = mul.cbegin(); arc != mul.cend(); ++arc) { 
              if (this->graph->is_mobility_edge(s, arc->second, arc->first)) { 
                mobiles.insert(*arc); 
              } else { 
                non_mobiles.insert(*arc);
              }
            } 
          };

          mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          std::set<edge_t> mobile_edges_x; 
          std::set<edge_t> non_mobile_edges_x;
          filter_lambda(x, mularcs_x, mobile_edges_x, non_mobile_edges_x);

          mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y);
          mularcs_y.erase(x);
          
          std::set<edge_t> mobile_edges_y; 
          std::set<edge_t> non_mobile_edges_y;
          filter_lambda(y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          if (non_mobile_edges_x.empty() && non_mobile_edges_y.empty()) {
          
            typedef std::pair<std::pair<vertex_t, structure::Mcolor>, size_t> colacr_t;
            utility::equivalence<colacr_t> equiv; 
            for (auto const & arc_x : mularcs_x) { 
              for (auto const & arc_y : mularcs_y) {
                mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
                if (color.size() > 0) {   
                  equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
                } 
              } 
            }  
            equiv.update();

            std::list<clone_t> clones; 
            std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.get_eclasses<std::set<colacr_t> >();   
            bool flag = false; 
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

              if (mularcs_left.size() == 1 && mularcs_right.size() == 1 && mularcs_left.begin()->second == mularcs_right.begin()->second) { 
                flag = true; 
                break;
              } else if ((mularcs_left.size() == 1) || (mularcs_right.size() == 1)) { 
                mularcs_t mularcs_mother; 
                mularcs_t mularcs_father;
              
                if (mularcs_left.size() == 1) { 
                  mularcs_mother = mularcs_left; 
                  mularcs_father = mularcs_right;
                } else if (mularcs_right.size() == 1) { 
                  mularcs_mother = mularcs_right; 
                  mularcs_father = mularcs_left;
                }

                bool sligshot = (mularcs_mother.size() == 1) && (mularcs_father.size() != 1) && this->graph->is_vec_T_consistent_color(mularcs_mother.cbegin()->second);
                for (auto arc = mularcs_father.cbegin(); arc != (mularcs_father.cend()) && sligshot; ++arc) {
                  sligshot = this->graph->is_vec_T_consistent_color(arc->second);
                }
                sligshot = sligshot && (mularcs_mother.cbegin()->second == mularcs_father.union_multicolors());   
                
                if (sligshot) { 
                  vertex_t const & mother = mularcs_mother.begin()->first;        

                  if (mother == Infty) { 
                    std::string pseudo_vertex = "o" + std::to_string(pseudo_infinity_vertex) + "o";
                    ++pseudo_infinity_vertex;
                    clone_t clone(arc_t(x, y), mularcs_father, edge_t(pseudo_vertex, mularcs_mother.cbegin()->second), true);  
                    clones.push_back(clone);
                  } else {
                    clone_t clone(arc_t(x, y), mularcs_father, *(mularcs_mother.cbegin()), false);  
                    clones.push_back(clone);
                  }  
                }
              } else { 
                assert(false);
              }
            } 

            if (!flag && is_good_clones(clones)) {
              for (auto const & clone : clones) {
                //std::cerr << "Do clone in first case" << std::endl;
                found = true;
                this->graph->apply(clone);
                ++number_rear;
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
*/

//experemental stage
template<class graph_t>
bool Algorithm<graph_t>::stage5_4() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
      number_rear = 0;
        utility::equivalence<vertex_t> classes = graph->split_on_components();
        std::map<vertex_t, std::set<vertex_t> > components = classes.get_eclasses<std::set<vertex_t> >();


      for (auto const & component : components) {
        bool flag = false; 

        std::set<vertex_t> const & verteces = component.second;
        std::unordered_set<vertex_t> processed;

        std::cerr << "Start process component" << std::endl;
        for (auto const v : verteces) {
          std::cerr << v << " "; 
        } 
        std::cerr << std::endl;

        std::map<mcolor_t, std::vector<arc_t> > how_many_edges;
        for (auto const & vertex : verteces) { 
          mularcs_t mularcs = graph->get_all_adjacent_multiedges_with_info(vertex);
          for (auto const & arc : mularcs) { 
            if (processed.count(arc.first) == 0 && graph->is_vec_T_consistent_color(arc.second)) {
              how_many_edges[arc.second].push_back(std::make_pair(vertex, arc.first));
            }             
          }
          processed.insert(vertex);
        }

        for (auto const & edges_with_color : how_many_edges) {
          if (edges_with_color.second.size() == 2) {
            auto const & first_edge = *edges_with_color.second.cbegin();
            auto const & second_edge = *(edges_with_color.second.cbegin() + 1);

            mularcs_t mularcs_x = graph->get_all_adjacent_multiedges_with_info(first_edge.first);
            mularcs_t mularcs_y = graph->get_all_adjacent_multiedges_with_info(first_edge.second);

            twobreak_t twobreak;
            //KILL ME FOR THIS CODE, but I write it very fast and with end of working day. 
            if (mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge, second_edge, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.first, first_edge.second, second_edge.second, second_edge.first, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.second, first_edge.first, second_edge.first, second_edge.second, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.second, first_edge.first, second_edge.second, second_edge.first, edges_with_color.first);
            } else if (mularcs_x.defined(second_edge.first) && !mularcs_x.defined(second_edge.second)
              && !mularcs_y.defined(second_edge.first) && mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge, second_edge, edges_with_color.first);
            } else if (!mularcs_x.defined(second_edge.first) && mularcs_x.defined(second_edge.second)
              && mularcs_y.defined(second_edge.first) && !mularcs_y.defined(second_edge.second)) {
              twobreak = twobreak_t(first_edge.first, first_edge.second, second_edge.second, second_edge.first, edges_with_color.first);
            }

            std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
            graph->apply(twobreak);
            ++number_rear;
            flag = true; 
            break;
          }
        }
      }

      if (number_rear != 0) { 
      isChanged = true;
    }
  } while (number_rear > 0); 

  return isChanged;   
}

typedef typename utility::equivalence<vertex_t> equiv_t;
typedef typename std::map<vertex_t, std::set<arc_t> > bridges_t;

template<class graph_t>
bool Algorithm<graph_t>::stage5_3() { 
  bool isChanged = false;
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0;
      // go over all T-consistent multicolors
    for(auto ic = graph->cbegin_vec_T_consistent_color(); ic != graph->cend_vec_T_consistent_color(); ++ic) {
      const auto& Q = *ic;
      bool repeat = true;
      while(repeat) {
        repeat = false;
        
        std::map<vertex_t, std::set<arc_t> > EC; // reg. edges between diff. connected components of color Q
        std::map<vertex_t, std::set<arc_t> > EI; // irreg. edges of color Q
        equiv_t CC = this->graph->split_on_components_with_color(*ic);
        get_specific_edges(EC, EI, CC, *ic); 
        
        std::unordered_set<vertex_t> processed;
        // reg. edges between diff. connected components of color QQ
        for (const auto &reg_edge : EC) {
          if (processed.count(reg_edge.first) != 0) { 
            continue;
          }

          if (EI[reg_edge.first].size() != 0 || (reg_edge.second.size() % 2 != 0)) {
            continue;
          }

          std::set<vertex_t> verteces;
          partgraph_t bridges;
          for (auto const & arcs : reg_edge.second) { 
            verteces.insert(arcs.second);
            bridges.insert(arcs.first, arcs.second);
          }

          bool found = true;
          std::unordered_set<vertex_t> marks;
          std::list<twobreak_t> breaks;
          //std::cerr << "NEW COMPONENT" << std::endl;
          for (auto const & arcs : reg_edge.second) {
            if (marks.count(arcs.second) == 0) {
              //std::cerr << "Start to check " << arcs.second << std::endl;
              marks.insert(arcs.second);

              mularcs_t mularcs = graph->get_all_adjacent_multiedges(arcs.second);
              
              vertex_t match = ""; 
              size_t count_match = 0;
              for (auto marcs = mularcs.cbegin(); marcs != mularcs.cend(); ++marcs) {
                if (verteces.count(marcs->first) != 0) { 
                  match = marcs->first;
                  ++count_match;
                }
              }

              if (count_match == 1) {
                count_match = 0;
                mularcs_t mularcs_v = graph->get_all_adjacent_multiedges(match);
                for (auto marcs = mularcs_v.cbegin(); marcs != mularcs_v.cend(); ++marcs) {
                  if (verteces.count(marcs->first) != 0) { 
                    ++count_match;
                  }
                }       

                if (count_match == 1) {
                  marks.insert(match);
                  assert(bridges.find(match) != bridges.end());
                  vertex_t u = bridges.find(match)->second;
                  //std::cerr << "perfect find " << u << " " << match << std::endl;
                  twobreak_t br2(arcs.first, arcs.second, u, match, Q);
                  breaks.push_back(br2);
                } else { 
                  found = false;
                }
              } else {
                found = false;
              }
            }

            if (!found) {
              break;
            }
          }

          if (found) {
            //std::cerr << "ALL IS GOOD, do 2break" << std::endl;
            for (auto const & br2 : breaks) {
              graph->apply(br2);
              ++number_rear;
            }
            repeat = true;
            break;
          }
          //std::cerr << "END COMPONENT" << std::endl;
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
void Algorithm<graph_t>::get_specific_edges(bridges_t & regular_edges, bridges_t & irregular_edges, 
    equiv_t & connected_components, mcolor_t const & color) const {
  for(vertex_t const &x : *this->graph) {
    mularcs_t const & mularcs = this->graph->get_all_adjacent_multiedges(x); 
    for(auto const & arc : mularcs) {    
      if (arc.second == color) { 
        if (arc.first == Infty) { 
          irregular_edges[connected_components[x]].insert(std::make_pair(x, arc.first)); 
        } else if (!connected_components.isequiv(x, arc.first)) {
          regular_edges[connected_components[x]].insert(std::make_pair(x, arc.first)); 
        } 
      } 
    }
  }
}
 
/*
template<class graph_t>
bool Algorithm<graph_t>::stage71() {
  size_t number_rear = 0; // number of rearrangements 

  for (auto const & x: *graph) {  
    mularcs_t const & mularcs = graph->get_all_adjacent_multiedges(x);
      
    bool found = false;
    for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
      vertex_t const & y = im->first; // Q == im->second - color of central edge
    
      if (y == Infty) {
        continue;
      }

      if (!graph->is_T_consistent_color(im->second)) { 
      mularcs_t mularcs_x = graph->get_all_adjacent_multiedges(x);
      mularcs_x.erase(y);

      mularcs_t mularcs_y = graph->get_all_adjacent_multiedges(y);
      mularcs_y.erase(x);
 
      size_t count = 0;
      bool is_y = false;
      std::pair<vertex_t, mcolor_t> pr;
      for (auto left = mularcs_x.begin(); left != mularcs_x.end(); ++left) {  
        if (!graph->is_T_consistent_color(left->second)) {
          pr = *left;
          ++count;
        } 
      }     

      for (auto right = mularcs_y.begin(); right != mularcs_y.end(); ++right) { 
        if (!graph->is_T_consistent_color(right->second)) {
          pr = *right;
          is_y = true;
          ++count;
        } 
      }

      if (count == 1 && pr.first != Infty) { 
        size_t count_f = graph->split_color(im->second).size(); 
        size_t count_s = graph->split_color(pr.second).size();

        if (count_f == count_s) { 
    //std::cerr << "Situation " << x << " " << y << " " << pr.first << std::endl;
          mularcs_t mul_x = graph->get_all_adjacent_multiedges_with_info(x, false);
          mularcs_t mul_y = graph->get_all_adjacent_multiedges_with_info(y, false);
          size_t count_central = calculate_cost(y, mul_x, mul_y);
          size_t count_another = 0;
          if (is_y) { 
            auto temp = graph->get_all_adjacent_multiedges_with_info(pr.first, false);
            temp.erase(y);
            count_another = calculate_cost(pr.first, graph->get_all_adjacent_multiedges_with_info(y, false), temp);
    } else {
            auto temp = graph->get_all_adjacent_multiedges_with_info(pr.first, false);
            temp.erase(y);
           count_another = calculate_cost(pr.first, graph->get_all_adjacent_multiedges_with_info(x, false), temp);
          } 

          if (count_central == std::min(count_central, count_another)) {
            mul_x.erase(y);
            std::vector<twobreak_t> history;  
            bool good = true;
            for (auto const &arc : mul_x) {
        vertex_t const & v = mularcs_y.get_vertex(arc.second); 
              if (!v.empty() && graph->is_vec_T_consistent_color(arc.second)) {
                  //std::cerr << "Two_break " << x << " " << arc.first << " " << y << " " << v << " " << genome_match::mcolor_to_name(arc.second) << std::endl;                
                  history.push_back(twobreak_t(x, arc.first, y, v, arc.second));
            //graph->apply(twobreak_t(x, arc.first, y, v, arc.second));
              } else { 
                good = false;
              }  
            }

            if (good) { 
              found = true;
              for (auto const & break2 : history) {
                graph->apply(break2);
                ++number_rear; 
              }
            } 
          } 
        }
      }
      }  
    } 
  }
 
  return (number_rear != 0);
}
*/

#endif