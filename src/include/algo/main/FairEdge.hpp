#ifndef PROCESS_FAIR_EDGE_HPP 
#define PROCESS_FAIR_EDGE_HPP 

namespace algo { 

template<class graph_pack_t>
struct ProcessFairEdge : public algo::AbsStage<graph_pack_t> {
  
  using mcolor_t = typename graph_pack_t::mcolor_type;
  
  using edge_t = typename graph_pack_t::edge_t;  
  using arc_t = typename graph_pack_t::arc_t; 
  using mularcs_t = typename graph_pack_t::mularcs_t; 
  using twobreak_t = typename graph_pack_t::twobreak_t;
  
  explicit ProcessFairEdge(size_t max_round)
  : algo::AbsStage<graph_pack_t>("Process fair edges", "fair_edge", max_round) 
  {
  }
  
  bool run(graph_pack_t & graph_pack) override;
  
private: 
  using set_arc_t = std::set<arc_t>;

  bool is_good_twobreaks(graph_pack_t const & graph, std::vector<twobreak_t> const & twobreaks) const;

private:
  DECL_LOGGER("ProcessFairEdge");
};

template<class graph_pack_t>
bool ProcessFairEdge<graph_pack_t>::run(graph_pack_t & graph_pack) { 
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

        if (y != Infty && graph_pack.is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!graph_pack.is_mobility_edge(x, y)) { 
          TRACE("See on non mobile edge " << x << " " << y)
          mularcs_t mularcs_x = graph_pack.get_all_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          mularcs_t mularcs_y; 
          if (y != Infty) {  
            mularcs_y = graph_pack.get_all_adjacent_multiedges_with_info(y);    
            mularcs_y.erase(x); 
          }

          /*SPLIT ALL EDGES ON MOBILE AND NON MOBILE*/
          set_arc_t mobile_edges_x; 
          set_arc_t non_mobile_edges_x;
          algo::split_by_mobile_property(graph_pack, x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          set_arc_t mobile_edges_y; 
          set_arc_t non_mobile_edges_y;
          algo::split_by_mobile_property(graph_pack, y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

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
            if (is_good_twobreaks(graph_pack, possible_twobreaks)) {
              TRACE("Start do two-breaks in case 1 (Create complete brackets)")
              for (auto const & twobreak : possible_twobreaks) {
                {
                  std::ostringstream out; 
                  out << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " 
                      << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " 
                      << cfg::get().mcolor_to_name(twobreak.get_mcolor());
                  TRACE(out.str())
                }
                found = true;
                graph_pack.apply(twobreak);
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
              std::vector<twobreak_t> included_twobreaks;

              for (auto twobreak = possible_twobreaks.cbegin(); (twobreak != possible_twobreaks.cend()) && !additional_color.empty(); ++twobreak) { 
                mcolor_t action_color = twobreak->get_mcolor();
                if (additional_color.includes(action_color) && graph_pack.multicolors.is_vec_T_consistent_color(action_color)) {
                  included_twobreaks.push_back(*twobreak);
                  additional_color = mcolor_t(additional_color, action_color, mcolor_t::Difference);
                } 
              }

              if (additional_color.empty() && is_good_twobreaks(graph_pack, included_twobreaks)) { 
                TRACE("Start do two-breaks in case 2 (Create min T-consistent color)")
                {
                  auto str_color = cfg::get().mcolor_to_name(additional_color); 
                  TRACE("Need to create " << str_color << " on edge")
                }

                for (auto const & twobreak : included_twobreaks) {
                  {
                    std::ostringstream out; 
                    out << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " 
                        << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " 
                        << cfg::get().mcolor_to_name(twobreak.get_mcolor());
                    TRACE(out.str())
                  }
                
                  graph_pack.apply(twobreak);
                  found = true;
                  ++number_rear;
                }

                if (found) { 
                  assert(graph_pack.get_all_multicolor_edge(x, y).empty() || graph_pack.multicolors.is_T_consistent_color(graph_pack.get_all_multicolor_edge(x, y)));
                } 
              } 
            } 
          }

          /*CASE 3: Process unique two-breaks*/          
          if (!found) {
            TRACE("Start do two-breaks in case 3 (unique two-breaks)")
    
            for (auto const & twobreak : possible_twobreaks) {
              size_t count = graph_pack.mobility_score(twobreak.get_arc(0), twobreak.get_mcolor(), twobreak.get_arc(1)) + 
                  graph_pack.mobility_score(twobreak.get_arc(1), twobreak.get_mcolor(), twobreak.get_arc(0));      

              auto scores = graph_pack.is_decrease_verteces_score(twobreak);
              if ((count == 0) && (scores.first > scores.second)) { 
                {
                  std::ostringstream out; 
                  out << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " 
                      << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " 
                      << cfg::get().mcolor_to_name(twobreak.get_mcolor());
                  TRACE(out.str())
                }
                graph_pack.apply(twobreak);
                found = true;
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

  return isChanged;
} 

template<class graph_pack_t>
bool ProcessFairEdge<graph_pack_t>::is_good_twobreaks(graph_pack_t const & graph_pack, 
  std::vector<twobreak_t> const & twobreaks) const {
  bool is_all_good = true;

  for (auto vec_color = graph_pack.multicolors.cbegin_vec_T_consistent_color(); vec_color != graph_pack.multicolors.cend_vec_T_consistent_color(); ++vec_color) {
    size_t count_diff = 0;
    mcolor_t vec_target_color = *vec_color;

    for (twobreak_t const & twobreak : twobreaks) { 
      mcolor_t local_color = twobreak.get_mcolor();
      if (vec_target_color.includes(local_color)) {
        vec_target_color = mcolor_t(vec_target_color, local_color, mcolor_t::Difference);
        ++count_diff;
      }
    }
      
    if (vec_target_color.empty() && (count_diff != 1)) is_all_good = false;
  }

  return is_all_good;
}

} 

#endif
