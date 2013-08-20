#ifndef STAGE3_1_H_
#define STAGE3_1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage3() {
  size_t number_indel_event = 0; 
  
  std::unordered_set<vertex_t > processed; 
  for (const auto &a1 : *graph) {  
    const vertex_t& a2 = graph->get_obverse_vertex(a1);
    Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(a1);

    if (graph->is_indel_vertex(a1) && (processed.count(a1) == 0) && graph->is_indel_vertex(a2) && mularcs.size() != 0)  {
      processed.insert(a1); 
      processed.insert(a2);

      Mcolor indel_color = mularcs.union_multicolors(); 
      Mcolor bar_indel_color = graph->get_complement_color(indel_color);
      std::set<Mcolor> split_indel = graph->split_color(indel_color, false);
      std::set<Mcolor> split_bar_indel = graph->split_color(bar_indel_color, false);
      assert(indel_color == graph->get_adjacent_multiedges(a2).union_multicolors());

      if (graph->is_vec_T_consistent_color(bar_indel_color) 
	|| (split_bad_colors && ((split_bar_indel.size() == 2) || (split_indel.size() == split_bar_indel.size())))) { 
	//std::cerr << "Insertion: " << a1 << " " << a2 << " color: " << genome_match::mcolor_to_name(bar_indel_color) << std::endl;
	for (const auto &col : split_bar_indel) {
	  graph->apply_ins_del(InsDel<Mcolor>(a1, a2, col, false));
	  insertions.insert(std::make_pair(std::make_pair(a1, a2), col));
	  ++number_indel_event;
	}
      } else if ((!graph->is_vec_T_consistent_color(bar_indel_color) && graph->is_vec_T_consistent_color(indel_color))
	|| (split_bad_colors && (split_indel.size() == 2) && (split_bar_indel.size() > 2))) { 
	//std::cerr << "Postponed deletion: " << a1 << " " << a2 << " color: " << genome_match::mcolor_to_name(bar_indel_color) << std::endl;
	graph->apply_ins_del(InsDel<Mcolor>(a1, a2, bar_indel_color, false), false);
	postponed_deletions.insert(std::make_pair(std::make_pair(a1, a2), bar_indel_color));
	++number_indel_event;
      } 
    } 
  }
  return (number_indel_event != 0); 
} 

/*template<class graph_t>
void Algorithm<graph_t>::remove_postponed_deletions() {
  for (auto edge = postponed_deletions.begin(); edge != postponed_deletions.end();) { 
    const vertex_t& a1 = edge->first.first; 
    const vertex_t& a2 = edge->first.second;

    Mcolor color = graph->get_adjacent_multiedges(a1).get_multicolor(a2);
    if (color == graph->get_complete_color()) { 
      //std::cerr << "Start worked with viewed edge " << a1 << " " << a2 << " it's complete we remove it" << std::endl;
      graph->apply_ins_del(InsDel<Mcolor>(a1, a2, edge->second, true), false);      
      std::set<Mcolor> colors = graph->split_color(graph->get_complement_color(edge->second), false);
      for (const auto &col : colors) {
	InsDel<Mcolor> good_deletion(a1, a2, col, true); 
	graph->apply_ins_del(good_deletion);
      }
      postponed_deletions.erase(edge++);
    } else {
      //std::cerr << a1 << " " << a2 << " " << genome_match::mcolor_to_name(edge->second) << std::endl; 
      ++edge;
    }  
  } 
}*/
#endif
