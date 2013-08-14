#ifndef STAGE3_1_H_
#define STAGE3_1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage3_1() {
  auto is_good_split = [=] (const std::set<Mcolor>& loc_colors) -> bool {
    for (const auto& c : loc_colors) {
	if (!graph->is_vec_T_consistent_color(c)) {return false;}
    } 
    return true;
  };

  size_t number_indel_event = 0; 
  
  std::unordered_set<vertex_t > processed; 
  for (const auto &a1 : *graph) {  
    const vertex_t& a2 = graph->get_obverse_vertex(a1);
    Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(a1);

    if (graph->is_indel_vertex(a1) && (processed.count(a1) == 0) && graph->is_indel_vertex(a2) && mularcs.size() != 0)  {
      //std::cerr << "Start worked with " << a1 << " " << a2;
      processed.insert(a1); 
      processed.insert(a2);

      Mcolor indel_color = mularcs.union_multicolors(); 
      Mcolor bar_indel_color = graph->get_complement_color(indel_color);
      std::set<Mcolor> split_indel = graph->split_color(indel_color, false);
      std::set<Mcolor> split_bar_indel = graph->split_color(bar_indel_color, false);
      assert(indel_color == graph->get_adjacent_multiedges(a2).union_multicolors());

      if (graph->is_vec_T_consistent_color(bar_indel_color) 
	|| (split_bad_colors && (split_bar_indel.size() == 2))) {// && (split_indel.size() != 2))) {
	//std::cerr << " past vec-TC-color. Done." << std::endl; 
	for (const auto &col : split_bar_indel) {
		InsDel<Mcolor> insertion(a1, a2, col, false);
		graph->apply_ins_del(insertion);
		++number_indel_event;
	}
      } else if ((!graph->is_vec_T_consistent_color(bar_indel_color) && graph->is_vec_T_consistent_color(indel_color))
	|| (split_bad_colors && (split_indel.size() == 2) && (split_bar_indel.size() != 2))) { 
	//std::cerr << " past TC-color. Add to viewed edges." << std::endl;
	InsDel<Mcolor> bad_insertion(a1, a2, bar_indel_color, false); 
	graph->apply_ins_del(bad_insertion, false);
	viewed_edges.push_back(bad_insertion);
	viewed_edges1.insert(std::make_tuple(a1, a2, bar_indel_color));
	viewed_edges2.insert(std::make_pair(a1, a2));	
	++number_indel_event;
      } else if (split_bad_colors && (split_indel.size() != 2) && (split_bar_indel.size() != 2)) { 
	InsDel<Mcolor> bad_insertion(a1, a2, bar_indel_color, false); 
	graph->apply_ins_del(bad_insertion, false);
	viewed_edges.push_back(bad_insertion);
	viewed_edges1.insert(std::make_tuple(a1, a2, bar_indel_color));
	viewed_edges2.insert(std::make_pair(a1, a2));
	++number_indel_event;
      }  
    } 
  }

  return (number_indel_event != 0); 
} 

template<class graph_t>
void Algorithm<graph_t>::remove_past_bad_colors() {
  for (auto it = viewed_edges.begin(); it != viewed_edges.end(); ++it) {
    const vertex_t& a1 = it->get_edge().first; 
    const vertex_t& a2 = it->get_edge().second;

    //std::cerr << "Start worked with viewed edge " << a1 << " " << a2;
    Mcolor color = graph->get_adjacent_multiedges(a1).get_multicolor(a2);
    if (color == graph->get_complete_color()) {

      graph->apply_ins_del(it->inverse(), false);      
      std::set<Mcolor> colors = graph->split_color(graph->get_complement_color(it->get_mcolor()), false);
      for (const auto &col : colors) {
	InsDel<Mcolor> good_deletion(a1, a2, col, true); 
	graph->apply_ins_del(good_deletion);
      }

      viewed_edges.erase(it++);
      --it;
      //std::cerr << " yes, it's complete we remove it" << std::endl;
    } else {//if (color.empty()) { 
       std::cerr << a1 << " " << a2 << " " << genome_match::mcolor_to_name(it->get_mcolor()) << " " << graph->is_T_consistent_color(it->get_mcolor()) << " " << graph->split_color(it->get_mcolor(), false).size() << std::endl;
    } /*else {*/
      //std::cerr << " NO, now is bad " << genome_match::mcolor_to_name(color) << std::endl;
    //}
  } 
}

#endif
