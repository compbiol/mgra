#ifndef STAGE2_H_ 
#define STAGE2_H_ 

template<class graph_t>
size_t Algorithm<graph_t>::is_mobility_edge_score(const vertex_t& x, const mcolor_t& color, const vertex_t& y) const {
  std::unordered_set<vertex_t> processed; 

  mularcs_t&& mularcs_x = graph->get_adjacent_multiedges_with_info(x);
  mularcs_x.erase(y);
  mularcs_t&& mularcs_y = graph->get_adjacent_multiedges_with_info(y);
  mularcs_y.erase(x);

  const auto canformQ_Lambda = [&] (const vertex_t& current, const mcolor_t& Q) -> size_t {
    processed.insert(current);
    if (current == Infty) {
      return canformQoo;
    }

    const mularcs_t& mularcs = graph->get_adjacent_multiedges_with_info(current); 
    size_t canform = 1;
    for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
      
      mcolor_t inter_color(Q, arc->second, mcolor_t::Intersection); 
      if (inter_color.size() > 0 && inter_color.size() < arc->second.size()) { 
        canform = 0;
      }
    }
    return canform;
  };

  return 0;
} 

template<class graph_t>
bool Algorithm<graph_t>::is_mobility_edge(const vertex_t& x, const mcolor_t& color, const vertex_t& y) const {
  bool mobilQ = false;

  if (postponed_deletions.defined(x, y)) {
    return false;
  } 

  mularcs_t mularcs_x;
  if (x != Infty) {
    mularcs_x = graph->get_adjacent_multiedges_with_info(x);
    mularcs_x.erase(y);
  } else {
    mobilQ = canformQoo;
  }  

  mularcs_t mularcs_y;
  if (y != Infty) {
    mularcs_y = graph->get_adjacent_multiedges_with_info(x);
    mularcs_y.erase(x);
  } else {
    mobilQ = canformQoo;
  } 

  const auto canformQ_Lambda = [&] (const vertex_t& current, const mcolor_t& Q) -> bool {
    if (current == Infty) {
      return canformQoo;
    }
      
    const mularcs_t& mularcs = graph->get_adjacent_multiedges_with_info(current); 
    bool canform = true;
    for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
      mcolor_t color(Q, arc->second, mcolor_t::Intersection); 
      if (color.size() > 0 && color.size() < arc->second.size()) { 
        canform = false;
      }
    }
    return canform;
  };
  
  for(auto ix = mularcs_x.cbegin(); (ix != mularcs_x.cend()) && !mobilQ; ++ix) { 
    mobilQ = canformQ_Lambda(ix->first, color);
  } 

  if (!mobilQ) { 
    for(auto iy = mularcs_y.cbegin(); (iy != mularcs_y.cend()) && !mobilQ; ++iy) { 
      mobilQ = canformQ_Lambda(iy->first, color);
    }
  } 

  return mobilQ;    
}


template<class graph_t>
bool Algorithm<graph_t>::is_mobility_edge(const vertex_t& x, const vertex_t& y) const {
  size_t count_non_mobile = 0;

  if (postponed_deletions.defined(x, y)) {
    return false;
  } 
  
  bool mobilQ = false;
  //can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
  //if return false, then Q cannot be formed
  //if true - who knows...
  const auto canformQ_Lambda = [&] (const vertex_t& current, const mcolor_t& Q) -> bool {
    if (current == Infty) {
      return canformQoo;
    }
    // color Q can be formed if some it adjacent multicolors form a partition of Q
    // OR 
    // if every intersection Q \cap QQ = \emptyset or QQ.
    const mularcs_t& mularcs = graph->get_adjacent_multiedges_with_info(current); 
    bool canform = true;
    for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
      mcolor_t color(Q, arc->second, mcolor_t::Intersection); 
      if (color.size() > 0 && color.size() < arc->second.size()) { 
        canform = false;
      }
    }
    return canform;
  };
    
  if (x == Infty) {
      mobilQ = canformQoo;
  } else {
    const mularcs_t& mularcs_x = graph->get_adjacent_multiedges_with_info(x);
    const auto& arcs = mularcs_x.equal_range(y); 
    for (auto jc = arcs.first; (jc != arcs.second) && !mobilQ; ++jc) { 
      if (graph->is_vec_T_consistent_color(jc->second)) { //cental sub-edge
        for(auto ix = mularcs_x.cbegin(); (ix != mularcs_x.cend()) && !mobilQ; ++ix) { 
          if (ix->first != y) { 
            mobilQ = canformQ_Lambda(ix->first, jc->second);
          } 
        } 
      } 
    } 
  } 

  if (!mobilQ) { 
    if (y == Infty) {
      mobilQ = canformQoo;
    } else {
      const mularcs_t& mularcs_y = graph->get_adjacent_multiedges_with_info(y); 
      const auto& arcs = mularcs_y.equal_range(x); 
      for (auto jc = arcs.first; (jc != arcs.second) && !mobilQ; ++jc) { 
        if (graph->is_vec_T_consistent_color(jc->second)) { //cental sub-edge
          for(auto ix = mularcs_y.cbegin(); (ix != mularcs_y.cend()) && !mobilQ; ++ix) { 
            if (ix->first != x) { 
              mobilQ = canformQ_Lambda(ix->first, jc->second);
            } 
          } 
        } 
      }
    } 
  }  

  return mobilQ; 
} 

template<class graph_t> 
bool Algorithm<graph_t>::stage22() { 
 bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
	
    for(const auto &x : *graph) {  
      const mularcs_t& mularcs = graph->get_adjacent_multiedges(x);

      if (graph->is_duplication_vertex(x)) {
        continue;
      }  

      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
	const vertex_t& y = im->first; // Q == im->second - color of central edge

        if (y == Infty || graph->is_duplication_vertex(y)) {
          continue;
        } 
 
	if (!is_mobility_edge(x, y)) {  // != 0) {
          const auto filter_Lambda = [&] (const vertex_t& s, const mularcs_t& mul, std::unordered_set<vertex_t>& mobiles, std::unordered_set<vertex_t>& non_mobiles) -> void {
            for (auto arc = mul.cbegin(); arc != mul.cend(); ++arc) { 
              if (this->is_mobility_edge(s, arc->second, arc->first)) { // == 0) { 
                mobiles.insert(arc->first); 
              } else { 
                non_mobiles.insert(arc->first);
              }
            } 
          };  

          const auto canformQ_Lambda = [&] (const vertex_t& x, const mcolor_t& Q) -> size_t {
            if (x == Infty) {
              return canformQoo;
            }
            const mularcs_t& mularcs = graph->get_adjacent_multiedges_with_info(x); 
            bool canform = true;
	    for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
              mcolor_t color(Q, arc->second, mcolor_t::Intersection); 
              if (color.size() > 0 && color.size() < arc->second.size()) { 
                canform = false;
              }
            }
            return canform;
          };


          mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(x);
	  mularcs_x.erase(y);

          std::unordered_set<vertex_t> mobile_vertex_x; 
          std::unordered_set<vertex_t> non_mobile_vertex_x;
          filter_Lambda(x, mularcs_x, mobile_vertex_x, non_mobile_vertex_x);
          
          mularcs_t mularcs_y; 
          if (y != Infty) {  
            mularcs_y = graph->get_adjacent_multiedges_with_info(y);	  
            mularcs_y.erase(x); 
          } 

          std::unordered_set<vertex_t> mobile_vertex_y; 
          std::unordered_set<vertex_t> non_mobile_vertex_y;
	  filter_Lambda(y, mularcs_y, mobile_vertex_y, non_mobile_vertex_y);

          if (non_mobile_vertex_x.empty() && non_mobile_vertex_y.empty()) {
            for (const auto &arc : mularcs_x) {
              if (graph->is_vec_T_consistent_color(arc.second)) {
                const vertex_t& v = mularcs_y.get_vertex(arc.second); 
                if (!v.empty()) {   
        	  graph->apply_two_break(twobreak_t(x, arc.first, y, v, arc.second));
	          found = true;
	          ++number_rear;
                } 
              }  
            } 
          } else if (y != Infty) {
            const auto get_min_addit_color = [&] (const mcolor_t& color) -> mcolor_t {
              mcolor_t min_color = graph->get_complete_color();
              for (auto col = graph->cbegin_Tconsistent_color(); col != graph->cend_Tconsistent_color(); ++col) {
                if (col->includes(color)) {
                  mcolor_t diff_color(*col, color, mcolor_t::Difference);
                  if (diff_color.size() < min_color.size()) {
                    min_color = diff_color;
                  } 
                } 
              } 
              return min_color;
            };
  

	    mcolor_t addit_color = get_min_addit_color(im->second);
            if (!addit_color.empty() && addit_color != graph->get_complete_color()) {
              std::vector<std::pair<arc_t, mcolor_t> > endpoints; 
              for (auto arc = mularcs_x.cbegin(); (arc != mularcs_x.cend()) && !addit_color.empty(); ++arc) { 
                if (mobile_vertex_x.count(arc->first) != 0) { 
                  vertex_t v = mularcs_y.get_vertex(arc->second); 
                  if (mobile_vertex_y.count(v) != 0 && arc->second.includes(addit_color)) {
                    endpoints.push_back(std::make_pair(std::make_pair(arc->first, v), arc->second));   
                    addit_color = mcolor_t(addit_color, arc->second, mcolor_t::Difference);
 	          }  
                }
              }   
              
              if (addit_color.empty()) { 
                for(const auto& endpoint : endpoints) {
                  graph->apply_two_break(twobreak_t(x, endpoint.first.first, y, endpoint.first.second, endpoint.second));
                  found = true;
                  ++number_rear; 
                }  
              }
 	    }  
           
	    if (y != Infty && !found && !graph->is_vec_T_consistent_color(im->second)) {
              const auto count_variant_Lambda = [&] (const arc_t& viewed, const mcolor_t& Q, const arc_t& remove) -> size_t { 
                size_t number_variant = 0; 
                if (viewed.first != Infty) {
                  mularcs_t&& end_f = graph->get_adjacent_multiedges_with_info(viewed.first);
                  end_f.erase(viewed.second);
	          end_f.erase(remove.first); 
	          end_f.erase(remove.second); 
	          for (auto arc = end_f.cbegin(); arc != end_f.cend(); ++arc) { 
		    if (canformQ_Lambda(arc->first, Q)) {
                      ++number_variant;
                    } 
                  } 
                } 

		if (viewed.second != Infty) {                
                  mularcs_t&& end_s = graph->get_adjacent_multiedges_with_info(viewed.second);
                  end_s.erase(viewed.first);
	          end_s.erase(remove.first); 
	          end_s.erase(remove.second);
	          for (auto arc = end_s.cbegin(); arc != end_s.cend(); ++arc) { 
                    if (canformQ_Lambda(arc->first, Q)) {
                      ++number_variant;
                    }
                  }   
                } 
 
                return number_variant; 
              };
          
              for (auto arc = mularcs_x.cbegin(); arc != mularcs_x.cend(); ++arc) { 
		const vertex_t& v = mularcs_y.get_vertex(arc->second); 
                if ((mobile_vertex_x.count(arc->first) != 0) && (mobile_vertex_y.count(v) != 0)) { 
		  size_t count = count_variant_Lambda(arc_t(x, arc->first), arc->second, arc_t(y, v));
                  count += count_variant_Lambda(arc_t(y, v), arc->second, arc_t(x, arc->first));  
		  if (count == 0) {
                    graph->apply_two_break(twobreak_t(x, arc->first, y, v, arc->second));
		    //std::cerr << " Do 2-break in third case " << x << " " << y << std::endl;
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
    } 
           
    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 
   
  return isChanged;
}

template<class graph_t> 
bool Algorithm<graph_t>::stage2() { 
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
	
    for(const auto &x : *graph) {  
      const mularcs_t& mularcs = graph->get_adjacent_multiedges(x);

      if (graph->is_duplication_vertex(x)) {
        continue;
      }  
      
      bool found = false;
      for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
	const vertex_t& y = im->first; // Q == im->second - color of central edge

	if (y != Infty && !graph->is_duplication_vertex(y) && !is_mobility_edge(x, y)) {
          mularcs_t&& mularcs_x = graph->get_adjacent_multiedges_with_info(x);
	  mularcs_x.erase(y);
 	  mularcs_t&& mularcs_y = graph->get_adjacent_multiedges_with_info(y);
	  mularcs_y.erase(x);
          bool is_perfect_fair_edge = true; 

#ifdef VERSION2
	  for (auto arc = mularcs_x.cbegin(); (arc != mularcs_x.cend()) && is_perfect_fair_edge; ++arc) { 
	    if (arc->first == Infty) {
              is_perfect_fair_edge = !canformQoo;
            } else if (!is_mobility_edge(x, arc->first)) { 
	      is_perfect_fair_edge = false; 
            }
          } 

          for (auto arc = mularcs_y.cbegin(); (arc != mularcs_y.cend()) && is_perfect_fair_edge; ++arc) { 
	    if (arc->first == Infty) {
              is_perfect_fair_edge = !canformQoo;
            } else if (!is_mobility_edge(y, arc->first)) { 
	      is_perfect_fair_edge = false; 
            }
          } 
#endif
 	  if (is_perfect_fair_edge) { 
            for (const auto &arc : mularcs_x) { 
              const vertex_t& v = mularcs_y.get_vertex(arc.second);
              if (!v.empty() && graph->is_vec_T_consistent_color(arc.second)) {   
#ifdef LOG_ENABLED
                std::cerr << " Sub-multiedge " << x << " " << y << std::endl;
                std::cerr << x << " " << arc.first << " " << y << " " << v << genome_match::mcolor_to_name(arc.second) << std::endl;
#endif
	        graph->apply_two_break(twobreak_t(x, arc.first, y, v, arc.second));
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
 
#endif
