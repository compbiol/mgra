#ifndef STAGE3_1_H_
#define STAGE3_1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage3_1() {
  size_t number_indel_event = 0; 
  size_t ins = 0;
  size_t del = 0;
  size_t not_good = 0;

  std::unordered_set<vertex_t > processed; 
  for (const auto &a1 : graph) {  
    const vertex_t& a2 = graph.get_obverse_vertex(a1);
    Mularcs<Mcolor> mul = graph.get_adjacent_multiedges(a1);

    if (graph.is_indel_vertex(a1) && (processed.count(a1) == 0) && graph.is_indel_vertex(a2) && mul.size() != 0)  {
      //std::cerr << "Start worked with " << a1 << " " << a2;
      processed.insert(a1); 
      processed.insert(a2);

      Mcolor indel_color = mul.union_multicolors(); 
      Mcolor bar_indel_color = graph.get_complement_color(indel_color);
      assert(indel_color == graph.get_adjacent_multiedges(a2).union_multicolors());

      if (graph.is_vec_T_color(bar_indel_color)) {
	//std::cerr << " past vec-TC-color. Done." << std::endl; 
	InsDel<Mcolor> insertion(a1, a2, bar_indel_color, false);
	graph.apply_ins_del(insertion);
	++number_indel_event;
	++ins; 
      } else {
	if (graph.is_vec_T_color(indel_color) || split_bad_colors) { 
	  //std::cerr << " past TC-color. Add to viewed edges." << std::endl;
	  InsDel<Mcolor> bad_insertion(a1, a2, bar_indel_color, false); 
	  graph.apply_ins_del(bad_insertion, false);
	  viewed_edges.push_back(bad_insertion);
	  ++del;
	  ++number_indel_event;
	}
      } 
    }
  }

  //std::cerr << "Attempt worked with " << processed.size() << " vertex" << std::endl;
  //std::cerr << "have insertion (insert vec-TC-color) " << ins << std::endl; 
  //std::cerr << "have deletion (insert TC-color) " << del << std::endl;

  return (number_indel_event != 0); 
} 

template<class graph_t>
void Algorithm<graph_t>::remove_past_bad_colors() {
  for (auto it = viewed_edges.begin(); it != viewed_edges.end(); ++it) {
    const vertex_t& a1 = it->get_edge().first; 
    const vertex_t& a2 = it->get_edge().second;

    //std::cerr << "Start worked with viewed edge " << a1 << " " << a2;
    Mcolor color = graph.get_adjacent_multiedges(a1).get_multicolor(a2);
    if (color == graph.get_complete_color()) {
      graph.apply_ins_del(it->inverse(), false);
      
      std::set<Mcolor> colors = graph.split_color(graph.get_complement_color(it->get_mcolor()));
      for (const auto &col : colors) {
	InsDel<Mcolor> good_deletion(a1, a2, col, true); 
	graph.apply_ins_del(good_deletion);
      }

      viewed_edges.erase(it++);
      --it;
      //++number_indel_event;				
      //std::cerr << " yes, it's complete we remove it" << std::endl;
    } //else {
      //std::cerr << " NO, now is bad " << genome_match::mcolor_to_name(color) << std::endl;
    //}
  } 
}

#endif
