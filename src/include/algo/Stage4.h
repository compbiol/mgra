#ifndef STAGE4_1_ 
#define STAGE4_1_ 

template<class graph_t>
bool Algorithm<graph_t>::stage4_td() { 
  size_t number_dupl = 0;  
  size_t vtc = 0; 
  size_t bvtc = 0; 
  size_t tc = 0; 

  std::unordered_set<vertex_t> processed;

/*  for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
	if (processed.find(*it) == processed.end()) {
		std::vector<vertex_t> duplication; 
		processed.insert(*it); 
		duplication.push_back(*it);

		vertex_t current = *it; 
		do {
			current = graph.get_obverse_vertex(*it);	
			if (graph.is_have_self_loop(*it)) {
				processed.insert(current);
				duplication.push_back(current);			
			} 
			if (graph.is_duplication_vertex(current)) { 
				processed.insert(current);
				duplication.push_back(current);
			} 

			if (graph.get_adjacent_multiedge(current).is_have_one_to_one_match()) {
				processed.insert(); 
				duplication.push_back(current);
			} 
		} while

		Mularcs<Mcolor> mularcs = graph.get_adjacent_multiedges(*it); 
		Mcolor color = mularcs.get_multicolor(*it); 
			
		if (graph.is_vec_T_color(color)) { 
			++vtc; 
		} else if (graph.is_vec_T_color(graph.get_complement_color(color))) {
			++bvtc; 
		} else { 
			++tc;
		}
	}	
  }*/

  /*for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
	if (processed.find(*it) == processed.end()) { 
		if (graph.is_duplication_vertex(*it)) { 
			if (graph.is_have_self_loop(*it)) {
				std::cerr << " Self loop vertex have ";
			}
			std::cerr << *it << " have Multidegree " << graph.get_adjacent_multiedges(*it).size() << std::endl;
			Mularcs<Mcolor> mularcs = graph.get_adjacent_multiedges(*it); 		
			for(auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) { 
				std::cerr << genome_match::mcolor_to_name(im->second) << " to " << im->first << std::endl;  
			}
		} 
	} 	
  } */ 
  return false; 
} 

template<class graph_t>
bool Algorithm<graph_t>::stage4_rtd() {
  size_t number_dupl = 0; 	
  //size_t vtc = 0; 
  size_t bvtc = 0; 
  size_t tc = 0; 

  std::unordered_set<vertex_t> processed;

  for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
	if (graph.is_have_self_loop(*it)) { 
		std::cerr << "vertex " << *it;
		std::vector<arc_t> duplication({std::make_pair(*it, *it)}); 
		
		Mularcs<Mcolor> mularcs = graph.get_adjacent_multiedges(*it); 
		Mcolor color = mularcs.get_multicolor(*it); 
		
		if (graph.is_vec_T_color(color)) { 
			std::cerr << " have self-loop of a vecTC"; 
			vertex_t current = graph.get_obverse_vertex(*it);
			bool flag = true; 
			while (flag) {
				flag = false; 
				Mularcs<Mcolor> current_mularcs = graph.get_adjacent_multiedges(current);
				for (auto im = current_mularcs.cbegin(); im != current_mularcs.cend(); ++im) { 
					if (im->second.how_much_includes(color) == 2) {
						flag = true; 
						duplication.push_back(std::make_pair(current, im->first));
						current = im->first; 
						break;
					}	
				} 
				if (flag) {
					current = graph.get_obverse_vertex(current);
				}
			} 		

			std::cerr << " go to " << current;
			
			Mularcs<Mcolor> current_mularcs = graph.get_adjacent_multiedges(current);

			std::unordered_set<vertex_t> count_included;
			for (auto im = current_mularcs.cbegin(); im != current_mularcs.cend(); ++im) { 
				if (im->second.includes(color)) {
		 			count_included.insert(im->first);
				}
			} 
				
			if (count_included.size() != 2) {
				std::cerr << " have " << count_included.size() << " choise in vertex " << std::endl;  
				continue;					
			}

			Mcolor first = current_mularcs.get_multicolor(*count_included.begin());			
			Mcolor second = current_mularcs.get_multicolor(*(++count_included.begin()));
			if (first == color && second != color) {
				duplication.push_back(std::make_pair(*count_included.begin(), current));
			} else if (first != color && second == color) {
				duplication.push_back(std::make_pair(*(++count_included.begin()), current));
			} else if (graph.is_vec_T_color(first) && !graph.is_vec_T_color(second)) {
				duplication.push_back(std::make_pair(*(++count_included.begin()), current));
			} else if (!graph.is_vec_T_color(first) && graph.is_vec_T_color(second)) {
				duplication.push_back(std::make_pair(*count_included.begin(), current));
			} else {
				std::cerr << " don't know what choose" << std::endl;
				continue;
			}

			TandemDuplication<Mcolor> dupl(duplication, color, true);
			graph.apply_tandem_duplication(dupl);
			++number_dupl;
			std::cerr << std::endl;
 		} else if (graph.is_vec_T_color(graph.get_complement_color(color))) {
			//vertex_t current = graph.get_obverse_vertex(*it);
			std::cerr << " have self-loop of a TC" << std::endl; 
			++bvtc;
		} else { 
			std::cerr << " all have self-loop of a TC" << std::endl; 
			++tc;
		}		
	} 
  }

  //std::cerr << "We have self-loops: " << vtc + bvtc + tc << std::endl; 
  //std::cerr << "vec-TC-color self-loops: " << vtc << std::endl; 
  std::cerr << "bar vec-TC-color self-loops: " << bvtc << std::endl; 
  std::cerr << "non vec-TC-color self-loops: " << tc << std::endl; 

  if (number_dupl != 0) {
  	return true; 
  } else {
	return false;
  }
} 

#endif
