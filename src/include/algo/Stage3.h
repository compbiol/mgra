#ifndef STAGE3_1_H_
#define STAGE3_1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage3_1() {
	size_t number_indel_event = 0; 
	size_t ins = 0;
	size_t del = 0;
	size_t not_good = 0;

	std::unordered_set<vertex_t > processed; 
	for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) { 
		const vertex_t& a1 = *it; 
		const vertex_t& a2 = graph.get_obverse_vertex(*it);

		if (!graph.is_indel_vertex(a1) || (processed.find(a1) != processed.end()) || !graph.is_indel_vertex(a2))  {
			continue; 
		}	

		//std::cerr << "Start worked with " << a1 << " " << a2;

		processed.insert(a1); 
		processed.insert(a2);

		// ADD CHECK ABOUT NOT HAVE EDGE TO OBVERSE VERTEX
		Mcolor indel_color = graph.get_adjacent_multiedges(a1).union_multicolors(); 
		Mcolor bar_indel_color = graph.get_complement_color(indel_color);
		assert(indel_color == graph.get_adjacent_multiedges(a2).union_multicolors());

		if (graph.is_vec_T_color(bar_indel_color)) {
			//std::cerr << " past vec-TC-color. Done." << std::endl; 
			InsDel<Mcolor> insertion(a1, a2, bar_indel_color, false);
			graph.apply_ins_del(insertion);
			++number_indel_event;
			++ins; 
		} else if (graph.is_vec_T_color(indel_color)) {
			//std::cerr << " past TC-color. Add to viewed edges." << std::endl;
			InsDel<Mcolor> bad_insertion(a1, a2, bar_indel_color, false); 
			graph.apply_ins_del(bad_insertion, false);
			viewed_edges.push_back(std::make_pair(bad_insertion, std::set<Mcolor>({indel_color})));
			++number_indel_event;
			++del;
		} else {
			//std::cerr << " have both complimentary TC-color. Not worked." << std::endl;
			++not_good;
		}  	
	}


	std::cerr << "Attempt worked with " << processed.size() << " vertex" << std::endl;
	std::cerr << "have insertion (insert vec-TC-color) " << ins << std::endl; 
	std::cerr << "have deletion (insert TC-color) " << del << std::endl;
	std::cerr << "not process (both TC-color) " << not_good << std::endl;

	if (number_indel_event != 0) {
		return true; 
	} else {
		return false;
	}
} 

template<class graph_t>
bool Algorithm<graph_t>::stage3_2() {
	size_t number_indel_event = 0; 

	for (auto it = viewed_edges.begin(); it != viewed_edges.end(); ++it) {
		const vertex_t& a1 = it->first.get_edge().first; 
		const vertex_t& a2 = it->first.get_edge().second;

		//std::cerr << "Start worked with viewed edge " << a1 << " " << a2;
		Mcolor color = graph.get_adjacent_multiedges(a1).get_multicolor(a2);
		if (color == graph.get_complete_color()) {
			graph.apply_ins_del(it->first.inverse(), false);
			for (auto im = it->second.cbegin(); im != it->second.cend(); ++im) {
				InsDel<Mcolor> good_deletion(a1, a2, *im, true); 
				graph.apply_ins_del(good_deletion);
			}
			viewed_edges.erase(it++);
			--it;
			++number_indel_event;				
			//std::cerr << " yes, it's complite we remove it" << std::endl;
		} else {
			//std::cerr << " NO, now is bad " << genome_match::mcolor_to_name(color) << std::endl;
		}
	} 

	if (number_indel_event != 0) {
		return true; 
	} else {
		return false;
	}
}

template<class graph_t>
bool Algorithm<graph_t>::stage6() { //less reliable
	size_t number_indel_event = 0; 
	std::unordered_set<vertex_t > processed; 

	for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) { 
		const vertex_t& a1 = *it; 
		const vertex_t& a2 = graph.get_obverse_vertex(*it);

		if (!graph.is_indel_vertex(a1) || (processed.find(a1) != processed.end()) || !graph.is_indel_vertex(a2))  {
			continue; 
		}	

		Mularcs<Mcolor> mul1 = graph.get_adjacent_multiedges(a1); 		
		Mularcs<Mcolor> mul2 = graph.get_adjacent_multiedges(a2); 	
		Mcolor indel_color = mul1.union_multicolors();
		Mcolor bar_indel_color = graph.get_complement_color(indel_color); 

		if (mul1.size() == 0 || (mul1.get_multicolors() != mul2.get_multicolors()) || graph.is_vec_T_color(indel_color) || graph.is_vec_T_color(bar_indel_color)) { 
			continue; 
		} 

		bool flag = true; 
		for (auto it = mul1.cbegin(); it != mul1.cend(); ++it) {
			if (!graph.is_vec_T_color(it->second)) {
				flag = false; 
			}
		}

		if (!flag) {
			continue;
		}

		processed.insert(a1); 
		processed.insert(a2);
	
		std::cerr << a1 << " " << a2 << " " << mul1.size() << " have all vec TC COLOR " << std::endl; 
		
		InsDel<Mcolor> bad_insertion(a1, a2, bar_indel_color, false); 
		graph.apply_ins_del(bad_insertion, false);
		viewed_edges.push_back(std::make_pair(bad_insertion, mul1.get_multicolors()));
		++number_indel_event;
	}

	if (number_indel_event != 0) {
		return true; 
	} else {
		return false;
	} 
} 
#endif
