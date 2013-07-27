#ifndef STAGE4_1_ 
#define STAGE4_1_ 

template<class graph_t>
bool Algorithm<graph_t>::stage4_td() { 
  size_t number_dupl = 0;  
  size_t vtc = 0; 
  size_t bvtc = 0; 
  size_t tc = 0;
  size_t not_tandem = 0; 

  std::unordered_set<vertex_t> processed;

  for (auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
    if (graph.is_duplication_vertex(*it) && !graph.is_have_self_loop(*it)) {
      std::cerr << "vertex " << *it;
      Mularcs<Mcolor> mularcs = graph.get_adjacent_multiedges(*it); 
      if (mularcs.find(graph.get_obverse_vertex(*it)) != mularcs.cend()) {
	std::cerr << " --> " << graph.get_obverse_vertex(*it); 
	Mcolor color = mularcs.get_multicolor(graph.get_obverse_vertex(*it));	
	if (graph.is_vec_T_color(color)) { 
	  std::cerr << " is a vecTC. Process" << std::endl; 
	  ++vtc; 
	  std::vector<arc_t> duplication({std::make_pair(*it, graph.get_obverse_vertex(*it))}); 
	  TandemDuplication<Mcolor> dupl(duplication, color, true, false);
	  graph.apply_tandem_duplication(dupl);	
	  ++number_dupl;
	} else if (graph.is_vec_T_color(graph.get_complement_color(color))) {
	  std::cerr << " is a bar vecTC. Not worked now" << std::endl;
	  ++bvtc; 
	} else { 
	  std::cerr << " both is a TC. Not worked now" << std::endl;
	  ++tc;
	}
      } else {
	//std::cerr << " don't worked now" << std::endl;	
	bool find = false; 
	for (auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) {
   	  if (im->first == Infty) {
	    continue;
	  }

	  if (im->second.is_one_to_one_match()) {
	    Mcolor color = im->second;	
	    vertex_t current = im->first;
	    std::vector<arc_t> duplication({std::make_pair(*it, current)}); 
	    
	    while (current != *it) {	
	      Mularcs<Mcolor> current_mularcs = graph.get_adjacent_multiedges(current);
	      bool flag = false; 
	      for (auto is = current_mularcs.cbegin(); is != current_mularcs.cend(); ++is) {
		if (is->first == Infty) {
		  continue;
		}

		if (is->second.how_much_includes(color) == 2) {
		  flag = true; 
		  duplication.push_back(std::make_pair(current, is->first));
		  current = is->first; 
		  break;
		}
	      }

	      
	      if (flag) {
		current = graph.get_obverse_vertex(current);
	      } else {
		break;
	      }
	    }
	    
	    if (current == *it) {
	      if (graph.is_vec_T_color(color)) { 
		std::cerr << "->(" << duplication.rbegin()->first << "," << duplication.rbegin()->second << ") is a vecTC. Process" << std::endl; 
		++vtc; 
		TandemDuplication<Mcolor> dupl(duplication, color, true, false);
		graph.apply_tandem_duplication(dupl);	
		++number_dupl;
	      } else if (graph.is_vec_T_color(graph.get_complement_color(color))) {
		std::cerr << "->(" << duplication.rbegin()->first << "," << duplication.rbegin()->second << ") is a bar vecTC. Not worked now" << std::endl;
		++bvtc; 
	      } else { 
		std::cerr << "->(" << duplication.rbegin()->first << "," << duplication.rbegin()->second << ") both is a TC. Not worked now" << std::endl;
		++tc;
	      }
	      find = true;	
	      break;
	    }
	  }
	}
	if (!find) {
	  std::cerr << " not tandem duplication" << std::endl;
	  ++not_tandem;
	}
      }	
    }
  } 
  std::cerr << "we have tandem duplication " << vtc + bvtc + tc << std::endl;
  std::cerr << "vec-TC-color tandem duplication " << vtc << std::endl;
  std::cerr << "TC-color tandem duplication " << bvtc << std::endl;
  std::cerr << "both TC-color tandem duplication " << tc << std::endl;
  std::cerr << "Not tandem duplication " << not_tandem << std::endl;

  if (number_dupl != 0) {
    return true; 
  } else {
    return false; 
  }
} 

template<class graph_t>
bool Algorithm<graph_t>::stage4_rtd() {
  size_t number_dupl = 0; 	
  size_t vtc = 0; 
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
	++vtc;
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

	std::cerr << " --> " << current;
			
	Mularcs<Mcolor> current_mularcs = graph.get_adjacent_multiedges(current);

	std::unordered_set<vertex_t> count_included;
	for (auto im = current_mularcs.cbegin(); im != current_mularcs.cend(); ++im) { 
	  if (im->second.includes(color)) {
	    count_included.insert(im->first);
	  }
	} 
	if (count_included.size() != 2 && 
	    (count_included.size() != 1 || color != current_mularcs.get_multicolor(*count_included.begin()))) {
	  std::cerr << " have " << count_included.size() << " choise in vertex " << std::endl;  
	  continue;					
	}

	if (count_included.size() == 1 && color == current_mularcs.get_multicolor(*count_included.begin())) {
	  duplication.push_back(std::make_pair(*count_included.begin(), current));
	  duplication.push_back(std::make_pair(*it, *count_included.begin()));	
	} else {	
	  Mcolor first = current_mularcs.get_multicolor(*count_included.begin());			
	  Mcolor second = current_mularcs.get_multicolor(*(++count_included.begin()));
	  if (first == color && second != color) {
	    duplication.push_back(std::make_pair(*count_included.begin(), current));
	    duplication.push_back(std::make_pair(*it, *count_included.begin()));
	  } else if (first != color && second == color) {
	    duplication.push_back(std::make_pair(*(++count_included.begin()), current));
	    duplication.push_back(std::make_pair(*it, *(++count_included.begin())));
	  } else if (graph.is_vec_T_color(first) && !graph.is_vec_T_color(second)) {
	    duplication.push_back(std::make_pair(*(++count_included.begin()), current));
	    duplication.push_back(std::make_pair(*it, *(++count_included.begin())));
	  } else if (!graph.is_vec_T_color(first) && graph.is_vec_T_color(second)) {
	    duplication.push_back(std::make_pair(*count_included.begin(), current));
	    duplication.push_back(std::make_pair(*it, *count_included.begin()));
	  } else {
	    std::cerr << " don't know what choose" << std::endl;
	    continue;
	  }
	}

	TandemDuplication<Mcolor> dupl(duplication, color, true, true);
	graph.apply_tandem_duplication(dupl);
	++number_dupl;
	std::cerr << " processed" << std::endl;
      } else if (graph.is_vec_T_color(graph.get_complement_color(color))) {
	std::cerr << " have self-loop bar vecTC. Not worked now" << std::endl; 
	++bvtc;
      } else { 
	std::cerr << " have self-loop both TC. Not worked now" << std::endl; 
	++tc;
      }		
    } 
  }

  std::cerr << "We have self-loops: " << vtc + bvtc + tc << std::endl; 
  std::cerr << "vec-TC-color self-loops: " << vtc << std::endl; 
  std::cerr << "bar vec-TC-color self loops: " << bvtc << std::endl; 
  std::cerr << "both is TC in self-loop: " << tc << std::endl; 

  if (number_dupl != 0) {
    return true; 
  } else {
    return false;
  }
} 

#endif
