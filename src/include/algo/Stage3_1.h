#ifndef STAGE3_1_H_
#define STAGE3_1_H_

template<class graph_t>
bool Algorithm<graph_t>::newstage3_1() {
	size_t number_indel_event = 0; 
	size_t ins = 0;
	size_t del = 0;

	std::cerr << "Start work with stage 3" << std::endl;

	std::unordered_set<vertex_t > processed; 
	for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) { 
		const vertex_t& a1 = *it; 
		const vertex_t& a2 = graph.get_obverse_vertex(*it);

		if (!graph.is_indel_vertex(a1) || (processed.find(a1) != processed.end()) || !graph.is_indel_vertex(a2))  {
			continue; 
		}	

		processed.insert(a1); 
		processed.insert(a2);

		Mcolor indel_color = graph.get_adjacent_multiedges(a1).union_multicolors(); 
		assert(indel_color == graph.get_adjacent_multiedges(a2).union_multicolors());
		Mcolor bar_indel_color = graph.get_complement_color(indel_color);

		if (graph.is_vec_T_color(bar_indel_color)) {
			++ins; 
			InsDel<Mcolor> insertion(a1, a2, bar_indel_color, false);
			graph.apply_ins_del(insertion);
			++number_indel_event;
		} else if (graph.is_vec_T_color(indel_color)) {
			++del;
			InsDel<Mcolor> bad_insertion(a1, a2, bar_indel_color, false); 
			graph.apply_ins_del(bad_insertion, false);
			viewed_edges.insert(std::make_pair(a1, bar_indel_color));
			++number_indel_event;
		}  	
	}


	std::cerr << "Worked with " << processed.size() << std::endl;
	std::cerr << "Have insertion " << ins << std::endl; 
	std::cerr << "Have deletion " << del << std::endl;

	if (number_indel_event != 0) {
		return true; 
	} else {
		return false;
	}
} 

template<class graph_t>
bool Algorithm<graph_t>::newstage3_2() {
	size_t number_indel_event = 0; 

	for (auto it = viewed_edges.cbegin(); it != viewed_edges.cend(); ++it) {
		const vertex_t& a1 = it->first; 
		const vertex_t& a2 = graph.get_obverse_vertex(it->first);

		// FIXME: IF WE CREATE DUPLICATION VERTEX? 
		Mcolor color = graph.get_adjacent_multiedges(a1).get_multicolor(a2);
		if (color == genome_match::get_complite_color()) {
			InsDel<Mcolor> bad_deletion(a1, a2, it->second, true); 
			graph.apply_ins_del(bad_deletion, false);
			assert(graph.is_vec_T_color(graph.get_complement_color(it->second)));
			InsDel<Mcolor> good_deletion(a1, a2, graph.get_complement_color(it->second), true); 
			graph.apply_ins_del(good_deletion);
			viewed_edges.erase(it);
			--it;
			++number_indel_event;				
		} 
	} 

	if (number_indel_event != 0) {
		return true; 
	} else {
		return false;
	}
} 
#endif
