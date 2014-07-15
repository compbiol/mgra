#ifndef STAGE3_H_
#define STAGE3_H_

template<class graph_t>
bool Algorithm<graph_t>::stage3() {
  size_t number_indel_event = 0; 
  
  std::unordered_set<vertex_t > processed; 
  for (auto const &a1 : *graph) {  
    vertex_t const & a2 = graph->get_obverse_vertex(a1);
    mularcs_t const & mularcs = graph->get_adjacent_multiedges(a1);

    if (graph->is_indel_vertex(a1) && (processed.count(a1) == 0) && graph->is_indel_vertex(a2) && mularcs.size() != 0)  {
      processed.insert({a1, a2}); 
      
      auto const & indel_color = mularcs.union_multicolors(); 
      auto const & bar_indel_color = graph->get_complement_color(indel_color);
      assert(indel_color == graph->get_adjacent_multiedges(a2).union_multicolors());
      graph->apply(insertion_t(a1, a2, bar_indel_color, false));
      ++number_indel_event;

      auto set_split_indel = graph->split_color(indel_color); 
      bool i_tc = false;
      size_t c_vtc_indel = 0; 
      for (auto const & col: set_split_indel) {  
        if (graph->is_vec_T_consistent_color(col)) {
          ++c_vtc_indel; 
        } else { 
          i_tc = true; 
        } 
      }  
      
      auto set_split_bar_indel = graph->split_color(bar_indel_color); 
      bool bi_tc = false;
      size_t c_vtc_bar_indel = 0;        
      for (auto const & col: set_split_bar_indel) {  
        if (graph->is_vec_T_consistent_color(col)) {
          ++c_vtc_bar_indel; 
        } else { 
          bi_tc = true; 
        }
      }  

      if (i_tc || c_vtc_bar_indel == std::min(c_vtc_indel, c_vtc_bar_indel)) { 
        insertions.insert(a1, a2);
      } else if (bi_tc || c_vtc_indel == std::min(c_vtc_indel, c_vtc_bar_indel)) { 
        graph->registrate_viewed_edge(a1, a2);
        postponed_deletions.insert(a1, a2);
      } else { 
        assert(false);
      } 
      
      /*size_t degree_split_bar_indel = graph->max_degree_split_color(bar_indel_color);
      size_t degree_split_indel = graph->max_degree_split_color(indel_color);
      if (degree_split_bar_indel == std::min(degree_split_bar_indel, degree_split_indel)) { 
        insertions.insert(a1, a2);
      } else if (degree_split_indel == std::min(degree_split_bar_indel, degree_split_indel)) { 
        graph->registrate_viewed_edge(a1, a2);
	postponed_deletions.insert(a1, a2);
      }*/
    } 
  }
  return (number_indel_event != 0); 
} 

template<class graph_t>
size_t Algorithm<graph_t>::check_postponed_deletions() const {
  size_t bad_postponed_deletion = 0;

  for (auto const & edge : postponed_deletions) {
    vertex_t const & a1 = edge.first; 
    vertex_t const & a2 = edge.second;

    mcolor_t const & color = graph->get_adjacent_multiedges(a1).get_multicolor(a2);
    if (color != graph->get_complete_color()) {    
      std::cerr << a1 << " " << a2 << std::endl;
      ++bad_postponed_deletion; 
    } 
  }

  return bad_postponed_deletion;
}
#endif
