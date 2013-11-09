#ifndef STAGE3_H_
#define STAGE3_H_

template<class graph_t>
bool Algorithm<graph_t>::stage3() {
  size_t number_indel_event = 0; 
  
  std::unordered_set<vertex_t > processed; 
  for (const auto &a1 : *graph) {  
    const vertex_t& a2 = graph->get_obverse_vertex(a1);
    const mularcs_t& mularcs = graph->get_adjacent_multiedges(a1);

    if (graph->is_indel_vertex(a1) && (processed.count(a1) == 0) && graph->is_indel_vertex(a2) && mularcs.size() != 0)  {
      processed.insert({a1, a2}); 
      
      const auto& indel_color = mularcs.union_multicolors(); 
      const auto& bar_indel_color = graph->get_complement_color(indel_color);
      assert(indel_color == graph->get_adjacent_multiedges(a2).union_multicolors());
      graph->apply_ins_del(insertion_t(a1, a2, bar_indel_color, false));
      ++number_indel_event;

#ifdef ROOT_LEAF
      if (!bar_indel_color.includes(graph->get_root_color())) { 
        insertions.insert(a1, a2);
      } else { 
        graph->registrate_viewed_edge(a1, a2);
	postponed_deletions.insert(a1, a2);
      }
#endif
 
      size_t degree_split_indel = graph->max_degree_split_color(indel_color);
      size_t degree_split_bar_indel = graph->max_degree_split_color(bar_indel_color); 	
      if (degree_split_bar_indel == std::min(degree_split_bar_indel, degree_split_indel)) { 
        insertions.insert(a1, a2);
      } else if (degree_split_indel == std::min(degree_split_bar_indel, degree_split_indel)) { 
        graph->registrate_viewed_edge(a1, a2);
	postponed_deletions.insert(a1, a2);
      }
    } 
  }
  return (number_indel_event != 0); 
} 

template<class graph_t>
size_t Algorithm<graph_t>::check_postponed_deletions() const {
  size_t bad_postponed_deletion = 0;

  for (const auto& edge : postponed_deletions) {
    const vertex_t& a1 = edge.first; 
    const vertex_t& a2 = edge.second;

    const auto& mularcs = graph->get_adjacent_multiedges(a1);
    const auto& pair = mularcs.equal_range(a2);
    mcolor_t color;
    for (auto it = pair.first; it != pair.second; ++it) {
      color = mcolor_t(color, it->second, mcolor_t::Union); 
    } 

    if (color != graph->get_complete_color()) {    
      std::cerr << a1 << " " << a2 << std::endl;
      ++bad_postponed_deletion; 
    } 
  }

  return bad_postponed_deletion;
}
#endif
