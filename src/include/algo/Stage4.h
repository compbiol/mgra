#ifndef STAGE4_H_ 
#define STAGE4_H_ 

template<class graph_t>
bool Algorithm<graph_t>::stage4_td() { 
  size_t number_dupl = 0;  
  std::unordered_set<vertex_t> processed;

  for (const auto &v : *graph) {
    if (graph->is_duplication_vertex(v) && !graph->is_have_self_loop(v)) {
      mularcs_t mularcs = graph->get_adjacent_multiedges(v); 
      bool find = false; 

      for (auto im = mularcs.cbegin(); im != mularcs.cend() && !find; ++im) {
	if (im->first != Infty && im->second.is_one_to_one_match()) {
	  vertex_t current = im->first;
	  auto const & color = im->second;	
	  std::vector<arc_t> duplication({std::make_pair(v, current)}); 
	 
	  bool is_go = true;    
	  while (v != graph->get_obverse_vertex(current) && is_go) {	
	    current = graph->get_obverse_vertex(current); 
	    mularcs_t current_mularcs = graph->get_adjacent_multiedges(current);
	    is_go = false; 
	    for (auto is = current_mularcs.cbegin(); (is != current_mularcs.cend()) && !is_go; ++is) {
	      if (is->first != Infty && is->second.how_much_includes(color) >= 2) {
      		is_go = true; 
      		duplication.push_back(std::make_pair(current, is->first));
      		current = is->first; 
	      }
	    }
 	  }

	  if (v == graph->get_obverse_vertex(current)) {
	    if (graph->is_vec_T_consistent_color(color)) { 
	      tandem_duplication_t dupl(duplication, color, true, false);
	      graph->apply(dupl);	
	      ++number_dupl;
	    } else if (!graph->is_vec_T_consistent_color(graph->get_complement_color(color)) && split_bad_colors) {
		auto colors = graph->split_color(color);
		for(const auto &col: colors) {
		  tandem_duplication_t dupl(duplication, col, true, false);
	          graph->apply(dupl);	
	          ++number_dupl;
		} 
	    } 
            find = true;	
          }
	}
      }
    }
  } 

  return (number_dupl != 0); 
} 

template<class graph_t>
bool Algorithm<graph_t>::stage4_rtd() {
  size_t number_dupl = 0; 	
  std::unordered_set<vertex_t> processed;

  for (const auto &v : *graph) {
    if (graph->is_have_self_loop(v)) { 
      //std::cerr << "vertex " << *it;
      std::vector<arc_t> duplication({std::make_pair(v, v)}); 
		
      mularcs_t mularcs = graph->get_adjacent_multiedges(v); 
      vertex_t current = graph->get_obverse_vertex(v);
      const auto& color = mularcs.get_multicolor(v); 
      bool flag = true; 

      while (flag) {
	flag = false; 
	mularcs_t current_mularcs = graph->get_adjacent_multiedges(current);
	for (auto im = current_mularcs.cbegin(); (im != current_mularcs.cend()) && !flag; ++im) { 
	  if (im->second.how_much_includes(color) == 2) {
	    flag = true; 
	    duplication.push_back(std::make_pair(current, im->first));
	    current = im->first; 
	  }	
	} 
	if (flag) {
	  current = graph->get_obverse_vertex(current);
	}
      } 		

      //std::cerr << " --> " << current;
			
      mularcs_t current_mularcs = graph->get_adjacent_multiedges(current);
      std::unordered_set<vertex_t> count_included;

      for (const auto &arc : current_mularcs) { 
	if (arc.second.includes(color)) {
	  count_included.insert(arc.first);
	}
      } 
      if (count_included.size() != 2 && 
	  (count_included.size() != 1 || color != current_mularcs.get_multicolor(*count_included.begin()))) {
	//std::cerr << " have " << count_included.size() << " choise in vertex " << std::endl;  
	continue;					
      }

      if (count_included.size() == 1 && color == current_mularcs.get_multicolor(*count_included.begin())) {
	duplication.push_back(std::make_pair(*count_included.begin(), current));
	duplication.push_back(std::make_pair(v, *count_included.begin()));	
      } else {	
	auto first = current_mularcs.get_multicolor(*count_included.begin());			
	auto second = current_mularcs.get_multicolor(*(++count_included.begin()));
	if (first == color && second != color) {
	  duplication.push_back(std::make_pair(*count_included.begin(), current));
	  duplication.push_back(std::make_pair(v, *count_included.begin()));
	} else if (first != color && second == color) {
	  duplication.push_back(std::make_pair(*(++count_included.begin()), current));
	  duplication.push_back(std::make_pair(v, *(++count_included.begin())));
	} else if (graph->is_vec_T_consistent_color(first) && !graph->is_vec_T_consistent_color(second)) {
	  duplication.push_back(std::make_pair(*(++count_included.begin()), current));
	  duplication.push_back(std::make_pair(v, *(++count_included.begin())));
	} else if (!graph->is_vec_T_consistent_color(first) && graph->is_vec_T_consistent_color(second)) {
	  duplication.push_back(std::make_pair(*count_included.begin(), current));
	  duplication.push_back(std::make_pair(v, *count_included.begin()));
	} else {
	  //std::cerr << " don't know what choose" << std::endl;
	  continue;
	}
      }

      if (graph->is_vec_T_consistent_color(color)) { 
	tandem_duplication_t dupl(duplication, color, true, true);
	graph->apply(dupl);
	++number_dupl;
	//std::cerr << " processed" << std::endl;
      } else if (!graph->is_vec_T_consistent_color(graph->get_complement_color(color)) && split_bad_colors) {
	auto colors = graph->split_color(color);
	for (const auto &col: colors) {
		tandem_duplication_t dupl(duplication, col, true, true);
		graph->apply(dupl);
		++number_dupl;
	} 
	//std::cerr << " processed" << std::endl;
      
      }		
    } 
  }

  return (number_dupl != 0); 
} 

template<class graph_t>
bool Algorithm<graph_t>::stage4_conv_to_td() { 
  size_t number_rear = 0;
  std::unordered_set<vertex_t> processed; 

  for (const auto &a1 : *graph) { 
    const vertex_t& a2 = graph->get_obverse_vertex(a1);
  
    if ((processed.count(a1) == 0) && graph->is_duplication_vertex(a1) && graph->is_duplication_vertex(a2)) {		
      processed.insert({a1, a2});
      mularcs_t mularcs_a1 = graph->get_adjacent_multiedges(a1); 
      mularcs_t mularcs_a2 = graph->get_adjacent_multiedges(a2);

      if (mularcs_a1.size() == 2 && mularcs_a1.size() == mularcs_a2.size()) { 
        mcolor_t color(mularcs_a1.cbegin()->second, mularcs_a1.crbegin()->second, mcolor_t::Intersection);
	if (((mularcs_a1.cbegin()->second == color && mularcs_a1.crbegin()->second != color)
	    || (mularcs_a1.cbegin()->second != color && mularcs_a1.crbegin()->second == color))
	    && ((mularcs_a2.cbegin()->second == color && mularcs_a2.crbegin()->second != color)
	    || (mularcs_a2.cbegin()->second != color && mularcs_a2.crbegin()->second == color))
           ) { 
	  const vertex_t& x = mularcs_a1.get_vertex(color); 
	  const vertex_t& y = mularcs_a2.get_vertex(color);
	  if (x != y && x != graph->get_obverse_vertex(y)) {
	  	if (graph->is_vec_T_consistent_color(color)) { 
			//std::cerr << a1 << " " << x << " " << a2 << " " << y << std::endl;
			twobreak_t t(a1, x, a2, y, color);
			graph->apply(t);
			++number_rear;
	  	} else if (!graph->is_vec_T_consistent_color(graph->get_complement_color(color)) && split_bad_colors) {
			//std::cerr << a1 << " " << x << " " << a2 << " " << y << std::endl;
			auto colors = graph->split_color(color);
			for (const auto& col: colors) { 
				twobreak_t t(a1, x, a2, y, col);
				graph->apply(t);
				++number_rear;
			}
		} 
	  } 
	}
      } 
    } 
  } 

  return (number_rear != 0);	
} 
#endif
