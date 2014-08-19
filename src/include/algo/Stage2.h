#ifndef STAGE2_H_ 
#define STAGE2_H_ 

template<class graph_t>
struct Algorithm<graph_t>::ProcessFairEdges : public Algorithm<graph_t>::Stage {
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  
  explicit ProcessFairEdges(std::shared_ptr<graph_t> const & graph, writer::Wdots<graph_t, ProblemInstance<mcolor_t> > w)
  : Stage(graph) 
  , writer(w)
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override {
    return "Process fair edges.";
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

  bool is_good_twobreaks(std::vector<twobreak_t> const & twobreaks) const;

private: //remove
writer::Wdots<graph_t, ProblemInstance<mcolor_t> > writer;
};

template<class graph_t>
bool Algorithm<graph_t>::ProcessFairEdges::do_action() { 
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

        if (y != Infty && this->graph->is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!this->graph->is_mobility_edge(x, y)) { 
          mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          mularcs_t mularcs_y; 
          if (y != Infty) {  
            mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y);    
            mularcs_y.erase(x); 
          }

          /*SPLIT ALL EDGES ON MOBILE AND NON MOBILE*/
          set_edge_t mobile_edges_x; 
          set_edge_t non_mobile_edges_x;
          split_by_mobile_property(x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          set_edge_t mobile_edges_y; 
          set_edge_t non_mobile_edges_y;
          split_by_mobile_property(y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          /*CREATE POSSIBLE MOBILE TWO-BREAKS*/
          size_t count_bad_edges = 0;
          std::vector<twobreak_t> possible_twobreaks;
          for (auto const &arc : mularcs_x) {
            vertex_t v = "";  
            if (y == Infty) {
              v = y; 
            } else { 
              v = mularcs_y.get_vertex(arc.second);
            }                  

            if (!v.empty() && (mobile_edges_x.count(arc) != 0) 
                  && (mobile_edges_y.count(std::make_pair(v, arc.second)) != 0)) { 
                possible_twobreaks.push_back(twobreak_t(x, arc.first, y, v, arc.second));
            } else { 
              ++count_bad_edges;
            }   
          }

          /*START TO WORK WITH MAIN ALGORITHM*/
          /*CASE 1: Create complete edge*/
          if (non_mobile_edges_x.empty() && non_mobile_edges_y.empty() && mobile_edges_x.size() == possible_twobreaks.size()) {
            assert(count_bad_edges == 0);
            if (is_good_twobreaks(possible_twobreaks)) {
              for (auto const & twobreak : possible_twobreaks) {
                //std::cerr << "Do two break in first case" << std::endl;
                //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                found = true;
                this->graph->apply(twobreak);
                ++number_rear;
              } 
            }
          } 

#ifdef PSEUDO_EDGE
          if (!found && non_mobile_edges_x.empty() && non_mobile_edges_y.empty()) {
            for (auto const & twobreak : possible_twobreaks) {
              size_t count = this->graph->mobility_score(twobreak.get_arc(0), twobreak.get_mcolor(), twobreak.get_arc(1)) + 
                  this->graph->mobility_score(twobreak.get_arc(1), twobreak.get_mcolor(), twobreak.get_arc(0));      
              if (count == 0) {
                found = true;
                this->graph->apply(twobreak);
                ++number_rear;
              } 
            }
          } 
#endif

          /*CASE 2: Create min T-consistent color*/          
          if (!found) {
            mcolor_t additional_color = get_min_addit_color_for_tc(this->graph->get_all_multicolor_edge(x, y));

            
            if (!additional_color.empty() && additional_color != this->graph->get_complete_color()) {
              std::vector<twobreak_t> included_twobreaks;

              for (auto twobreak = possible_twobreaks.cbegin(); (twobreak != possible_twobreaks.cend()) && !additional_color.empty(); ++twobreak) { 
                mcolor_t action_color = twobreak->get_mcolor();
                if (additional_color.includes(action_color) && this->graph->is_vec_T_consistent_color(action_color)) {
                  included_twobreaks.push_back(*twobreak);
                  additional_color = mcolor_t(additional_color, action_color, mcolor_t::Difference);
                } 
              }

              if (additional_color.empty() && is_good_twobreaks(included_twobreaks)) { 
                for (auto const & twobreak : included_twobreaks) {
                  //std::cerr << " Do 2-break in second case " << std::endl;
                  //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                  this->graph->apply(twobreak);
                  found = true;
                  ++number_rear;
                }
              } 
            } 
          } 

          /*CASE 3: Process unique two-breaks*/          
          if (!found) {
            for (auto const & twobreak : possible_twobreaks) {
              size_t count = this->graph->mobility_score(twobreak.get_arc(0), twobreak.get_mcolor(), twobreak.get_arc(1)) + 
                  this->graph->mobility_score(twobreak.get_arc(1), twobreak.get_mcolor(), twobreak.get_arc(0));      

              auto scores = this->graph->is_decrease_verteces_score(twobreak);

#ifdef PSEUDO_EDGE
              if ((count == 0) && (scores.first > scores.second)) {
#else
              if ((count == 0) && (scores.first > scores.second)) { 
#endif
                //std::cerr << " Do 2-break in third case " << std::endl;
                //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                this->graph->apply(twobreak);
                found = true;
                ++number_rear;
                   
              } 
            }
          }

          /*CASE 4: Create edge with tip*/          
          /*
          if (!found && !non_mobile_edges_x.empty() && !non_mobile_edges_y.empty()) {
            for (auto edge = non_mobile_edges_x.begin(); edge != non_mobile_edges_x.end() && !found; ++edge) { 
              if (this->graph->degree_vertex(edge->first) == 2) {
                mularcs_t mularcs_p = this->graph->get_all_adjacent_multiedges(edge->first); 
                mularcs_p.erase(x);

                if (this->graph->is_vec_T_consistent_color(mularcs_p.begin()->second) && 
                  non_mobile_edges_y.count(edge_t(mularcs_p.begin()->first, edge->second)) != 0) { 
                  mcolor_t target_color = mularcs_p.begin()->second;

                  //std::cerr << "central edge " << x << " " << y << std::endl;
                  //std::cerr << "Edge (p, q): " << edge->first << " " << mularcs_p.begin()->first << " " << genome_match::mcolor_to_name(target_color) << std::endl;
                  //std::cerr << "Edge (x, p) and (y, q) have color: " << genome_match::mcolor_to_name(edge->second) << std::endl;

                  for (auto const & twobreak : possible_twobreaks) { 
                    mcolor_t action_color = twobreak.get_mcolor(); 
                    if (target_color.includes(action_color)) { 
                      //std::cerr << " Do 2-break in fourth case " << std::endl;
                      std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                      this->graph->apply(twobreak);
                      ++number_rear;
                      found = true; 
                    } 
                  }
                }
              }
            } 
          } */

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
bool Algorithm<graph_t>::ProcessFairEdges::is_good_twobreaks(std::vector<twobreak_t> const & twobreaks) const {
  bool is_all_good = true;

  for (auto vec_color = this->graph->cbegin_vec_T_consistent_color(); vec_color != this->graph->cend_vec_T_consistent_color(); ++vec_color) {
    size_t count_diff = 0;
    mcolor_t vec_target_color = *vec_color;

    for (twobreak_t const & twobreak : twobreaks) { 
      mcolor_t local_color = twobreak.get_mcolor();
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
