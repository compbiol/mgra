#ifndef STAGE2_H_ 
#define STAGE2_H_ 

template<class graph_t>
struct Algorithm<graph_t>::ProcessFairEdges : public Algorithm<graph_t>::Stage {
  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  
  explicit ProcessFairEdges(std::shared_ptr<graph_t> const & graph)
  : Stage(graph) 
  {
  }
  
  bool do_action() override;
  
  std::string get_name() override {
    return "Process fair edges.";
  }

private: 
  size_t create_complete_edge();
  size_t create_nearest_T_consitent_color();
  size_t process_edge_with_unique_mobility_order();
  bool is_good_twobreaks(std::list<twobreak_t> const & twobreaks) const;
};

template<class graph_t>
bool Algorithm<graph_t>::ProcessFairEdges::do_action() { 
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
	
    for(vertex_t const & x : *this->graph) {  
      mularcs_t const & mularcs = this->graph->get_adjacent_multiedges(x);

      if (this->graph->is_duplication_vertex(x)) {
        continue;
      }  

      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
        vertex_t const & y = im->first; // Q == im->second - color of central edge

        if (y != Infty && this->graph->is_duplication_vertex(y)) {
          continue;
        } 
 
        if (!this->graph->is_mobility_edge(x, y) && im->second != this->graph->get_complete_color()) { 
          //std::cerr << "Work with non-mobile edge " << x << " " << y << std::endl;

          auto const filter_lambda = [&] (vertex_t const & s, mularcs_t const & mul, std::set<edge_t>& mobiles, std::set<edge_t>& non_mobiles) -> void {
            for (auto arc = mul.cbegin(); arc != mul.cend(); ++arc) { 
              if (this->graph->is_mobility_edge(s, arc->second, arc->first)) { 
                mobiles.insert(*arc); 
              } else { 
                non_mobiles.insert(*arc);
              }
            } 
          };

          mularcs_t mularcs_x = this->graph->get_adjacent_multiedges_with_info(x);
          mularcs_x.erase(y);

          std::set<edge_t> mobile_edges_x; 
          std::set<edge_t> non_mobile_edges_x;
          filter_lambda(x, mularcs_x, mobile_edges_x, non_mobile_edges_x);
          
          mularcs_t mularcs_y; 
          if (y != Infty) {  
            mularcs_y = this->graph->get_adjacent_multiedges_with_info(y);	  
            mularcs_y.erase(x); 
          } 

          std::set<edge_t> mobile_edges_y; 
          std::set<edge_t> non_mobile_edges_y;
          filter_lambda(y, mularcs_y, mobile_edges_y, non_mobile_edges_y);

          if (non_mobile_edges_x.empty() && non_mobile_edges_y.empty()) {
            size_t count_bad_edges = 0;
            std::list<twobreak_t> twobreaks; 

            for (auto const &arc : mularcs_x) {
              if (this->graph->is_vec_T_consistent_color(arc.second)) {
                vertex_t v = mularcs_y.get_vertex(arc.second); 
                if (v.empty() && y != Infty) { 
                  ++count_bad_edges;
                } else { 
                  if (y == Infty) {
                    v = y; 
                  }                  
                  twobreaks.push_back(twobreak_t(x, arc.first, y, v, arc.second));
                } 
              } else { 
                ++count_bad_edges;
              }    
            } 

            if (count_bad_edges == 0 && is_good_twobreaks(twobreaks)) {
              for (auto const & twobreak : twobreaks) {
                //std::cerr << "Do two break in first case" << std::endl;
                //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                found = true;
                this->graph->apply(twobreak);
                ++number_rear;
              } 
            }
          } 

          if (!found) { 
            auto const get_min_addit_color = [&] (mcolor_t const & color) -> mcolor_t {
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
            };
  

            std::list<twobreak_t> twobreaks; 
            mcolor_t addit_color = get_min_addit_color(im->second);
            
            if (!addit_color.empty() && addit_color != this->graph->get_complete_color()) {
              for (auto arc = mularcs_x.cbegin(); (arc != mularcs_x.cend()) && !addit_color.empty(); ++arc) { 
                if (mobile_edges_x.count(*arc) != 0 && this->graph->is_vec_T_consistent_color(arc->second)) { 
                  vertex_t v; 
                  if (y != Infty) {
                    v = mularcs_y.get_vertex(arc->second); 
                  } else {
                    v = y; 
                  } 
                  if ((y == Infty || mobile_edges_y.count(edge_t(v, arc->second)) != 0) && addit_color.includes(arc->second) && this->graph->is_vec_T_consistent_color(arc->second)) {
                    twobreaks.push_back(twobreak_t(x, arc->first, y, v, arc->second));
                    addit_color = mcolor_t(addit_color, arc->second, mcolor_t::Difference);
 	                } 
                }
              }   
              
              if (addit_color.empty() && is_good_twobreaks(twobreaks)) { 
                for (twobreak_t const & twobreak : twobreaks) {
                  //std::cerr << "Do two break in second case" << std::endl;
                  //std::cerr << twobreak.get_vertex(0) << " " << twobreak.get_vertex(1) << " " << twobreak.get_vertex(2) << " " << twobreak.get_vertex(3) << " " << genome_match::mcolor_to_name(twobreak.get_mcolor()) << std::endl;
                  this->graph->apply(twobreak);
                  found = true;
                  ++number_rear; 
                }  
              }
 	          }  
          }

          if (!found && !this->graph->is_vec_T_consistent_color(im->second)) {
            auto const count_variant_Lambda = [&] (arc_t const & viewed, mcolor_t const & Q, arc_t const & remove) -> size_t { 
              size_t number_variant = 0; 
              if (viewed.first != Infty) {
                mularcs_t&& end_f = this->graph->get_adjacent_multiedges_with_info(viewed.first);
                end_f.erase(viewed.second);
                end_f.erase(remove.first); 
                end_f.erase(remove.second); 
                for (auto arc = end_f.cbegin(); arc != end_f.cend(); ++arc) { 
                  if (this->graph->canformQ(arc->first, Q)) {
                    ++number_variant;
                  } 
                } 
              } 

	            if (viewed.second != Infty) {                
                mularcs_t&& end_s = this->graph->get_adjacent_multiedges_with_info(viewed.second);
                end_s.erase(viewed.first);
                end_s.erase(remove.first); 
                end_s.erase(remove.second);
                for (auto arc = end_s.cbegin(); arc != end_s.cend(); ++arc) { 
                  if (this->graph->canformQ(arc->first, Q)) {
                    ++number_variant;
                  }
                }   
              } 

              return number_variant; 
            };
        
            for (auto arc = mularcs_x.cbegin(); arc != mularcs_x.cend(); ++arc) { 
              vertex_t v; 
              if (y != Infty) {
                v = mularcs_y.get_vertex(arc->second); 
              } else { 
                v = y;
              } 
              
              if ((mobile_edges_x.count(*arc) != 0) && (mobile_edges_y.count(edge_t(v, arc->second)) != 0)) { 
	              size_t count = count_variant_Lambda(arc_t(x, arc->first), arc->second, arc_t(y, v));
                if (y != Infty) {
                  count += count_variant_Lambda(arc_t(y, v), arc->second, arc_t(x, arc->first));    
                }
	              if (count == 0 && this->graph->is_vec_T_consistent_color(arc->second)) {
                  //std::cerr << " Do 2-break in third case " << std::endl;
                  this->graph->apply(twobreak_t(x, arc->first, y, v, arc->second));
                  /*if (x == "261t" && arc->first == "1241h" && y == "702h" && v == "1158h") { 
                    std::cerr << genome_match::mcolor_to_name(this->graph->get_edge_multicolor(x, y)) 
                    << " " << genome_match::mcolor_to_name(this->graph->get_edge_multicolor(y, v)) << std::endl; 
                    exit(1);
                  }*/
	                //std::cerr << x << " " << arc->first << " " << y << " " << v << " " << genome_match::mcolor_to_name(arc->second) << std::endl;                
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
   
  return isChanged;
}

template<class graph_t>
size_t Algorithm<graph_t>::ProcessFairEdges::create_complete_edge() {
  size_t number_action = 0; 
  return number_action;
}

template<class graph_t>
size_t Algorithm<graph_t>::ProcessFairEdges::create_nearest_T_consitent_color() {
  size_t number_action = 0; 
  return number_action;
}

template<class graph_t>
size_t Algorithm<graph_t>::ProcessFairEdges::process_edge_with_unique_mobility_order() {
  size_t number_action = 0; 
  return number_action;
}
  
template<class graph_t>
bool Algorithm<graph_t>::ProcessFairEdges::is_good_twobreaks(std::list<twobreak_t> const & twobreaks) const {
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
