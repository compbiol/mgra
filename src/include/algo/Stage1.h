#ifndef STAGE1_H_
#define STAGE1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage1() {
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
    for (const auto& v: graph) {  
      if (graph.is_simple_vertex(v)) { 
	path_t path({v});
	std::unordered_set<vertex_t> processed({v, Infty}); // we count oo as already processed
	Mularcs<Mcolor> current = graph.get_adjacent_multiedges(v);

	for(auto im = current.cbegin(); im != current.cend(); ++im) {	
	  bool is_next = (im == current.cbegin()); 
	  vertex_t current = find_simple_path(path, processed, v, im->first, is_next);
	  if (current == v) { 
	    break; // got a cycle from x to x, cannot extend it 
	  }  		    
	}

	//if (!((*path.begin() == *path.rbegin()) && (*path.begin() != Infty) 
	//    && (graph.is_duplication_vertex(*path.begin()) || graph.is_indel_vertex(*path.begin())) 
	//      && (path.size() % 2 == 0))) {
	  number_rear += process_simple_path(path);
	//}
      }
    } 

    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 
   
  return isChanged;
} 

template<class graph_t>
vertex_t Algorithm<graph_t>::find_simple_path(path_t& path, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, bool is_next) { 
  vertex_t previous  = prev;
  vertex_t current = cur;
  bool stop = true; 

  while (stop) {
    stop = false; 

    if (is_next) { 
      path.push_front(current);
    } else { 
      path.push_back(current);
    } 

    if (processed.find(current) == processed.end() && !graph.is_duplication_vertex(current) && !graph.is_indel_vertex(current)) {     
      processed.insert(current);
      Mularcs<Mcolor> new_edges = graph.get_adjacent_multiedges(current);
      Mcolor previous_color = new_edges.get_multicolor(previous); 
      new_edges.erase(previous);
    
      if (new_edges.size() == 1 && graph.get_complement_color(previous_color) == new_edges.cbegin()->second) {
	  if (split_bad_colors) {
	    previous = current;
	    current = new_edges.cbegin()->first;  
	    stop = true;
	  } else {	
	    if (graph.is_T_consistent_color(new_edges.cbegin()->second) && graph.is_T_consistent_color(previous_color)) { 
	      previous = current;
	      current = new_edges.cbegin()->first; 
	      stop = true;
	    } 
	  }  
      } 
    }
  }

  return current;
} 

template<class graph_t>
size_t Algorithm<graph_t>::process_simple_path(path_t& path) {
  size_t number_rear = 0;

  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
    /*std::cerr << std::endl << "Processing a path of length " << path.size() - 1 << std::endl;
      std::cerr << "path:\t" << *path.begin();
      for(auto ip = ++path.begin(); ip != path.end(); ++ip) {
      std::cerr << " -- " << *ip;
      }
      std::cerr << std::endl;*/

    Mcolor process_color; 
    if (split_bad_colors) {
      Mcolor first = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
      Mcolor second = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*(++++path.begin()));
      if (graph.split_color(first).size() <= graph.split_color(second).size()) {
	process_color = first;
      } else { 
	process_color = second;
      }
    } else { 
      Mcolor first = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
      Mcolor second = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*(++++path.begin()));
      if (graph.is_vec_T_color(first)) {
	process_color = first;
      } else { 
	process_color = second;
      }
    }  

    if ((path.size() % 2 != 0) && (*path.begin() != *path.rbegin())) {
      //std::cerr << "... ";

      if (process_color != graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin())) { 
	path.erase(path.begin());
	//std::cerr << "left";
      } else {
	path.erase(--path.end());
	//std::cerr << "right";
      }
      //std::cerr << " end removed" << std::endl;
    }

    /*if (*path.begin() == Infty && *path.rbegin() == Infty ) {
      std::cerr << "... affecting two chromosome ends" << std::endl;
      } else if( *path.begin() == Infty || *path.rbegin() == Infty ) {
      std::cerr << "... affecting a chromosome end" << std::endl;
      }*/

    if (*path.begin() == *path.rbegin()) {
      if (path.size() % 2 == 0) {
	//assert(*path.begin() == Infty);
	//if (*path.begin() != Infty) {
	//std::cerr << "ERROR: Semi-cycle w/o infinity! " << *path.begin() << std::endl;
	//exit(1);
	//}

	if (process_color == graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin())) { 
	  //std::cerr << "... semi-cycle, fusion applied" << std::endl;
	  const vertex_t& self_v = *(path.begin());
	  const vertex_t& x0 = *(++path.begin());
	  const vertex_t& y0 = *(++path.rbegin());

	  Mularcs<Mcolor> mul = graph.get_adjacent_multiedges(x0, split_bad_colors);
	  auto colors = mul.equal_range(self_v);
	  for (auto it = colors.first; it != colors.second; ++it) { 
	    graph.apply_two_break(TwoBreak<Mcolor>(self_v, x0, self_v, y0, it->second));
	    ++number_rear;
	  } 
	
	  path.erase(--path.end());
	  *path.begin() = y0;
	} else {
	  //std::cerr << "... semi-cycle, fission applied" << std::endl;
	  const vertex_t& self_v = *(path.begin());
	  const vertex_t& y0 = *(++path.rbegin());
	  const vertex_t& y1 = *(++++path.rbegin());

	  Mularcs<Mcolor> mul = graph.get_adjacent_multiedges(y0, split_bad_colors);
	  auto pair = mul.equal_range(y1);
	  for (auto it = pair.first; it != pair.second; ++it) { 
 	    graph.apply_two_break(TwoBreak<Mcolor>(y0, y1, self_v, self_v, it->second));
	    ++number_rear;
	  }
        
	  path.erase(--path.end());
	  *path.rbegin() = self_v;
	}

	if (path.size() < 4) { 
	  return number_rear;
	}
      } //else { 
	//std::cerr << "... cycle" << std::endl;
      //} 
    }

    Mcolor Q = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
    while (process_color != Q) {
      // multicolor of (z1,z2). N.B.: x2 is NOT oo
      //std::cerr << "... multicolors of first and second multiedges: ";
      if (*path.begin() == *path.rbegin()) {
	//std::cerr << "... rotating" << std::endl;
	path.push_back(*path.begin());
	path.erase(path.begin());
      } else {
	if (*path.begin() == Infty && *path.rbegin() != Infty) {
	  //std::cerr << "... flipping" << std::endl;
	  for(auto ip = ++path.begin();ip != path.end();) {
	    path.push_front(*ip);
	    path.erase(ip++);
	  }
	}
	if (*path.rbegin() == Infty) {
	  //std::cerr << "... extending beyond oo" << std::endl;
	  path.push_back(Infty);
	  path.erase(path.begin());
	} else {
	  //std::cerr << "... truncating ??" << std::endl;
	  path.erase(path.begin());
	  path.erase(--path.end());
	  if (path.size() < 4) { 
	    return number_rear;
	  }
	}
      }
      Q = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
    }

    // x1 -- x2 -- x3 -- ... -- x2k
    // results in:
    // x2 == x3   x4 == x5 ... x(2k-2) == x(2k-1)   x1 -- x2k

    auto z3 = path.begin();
    auto z0 = z3++;
    auto z1 = z3++;
    auto z2 = z3++;
    std::set<Mcolor> colors = graph.split_color(Q);

    while(z3 != path.end()) {
      for (const auto &col: colors) {
	graph.apply_two_break(TwoBreak<Mcolor>(*z0, *z1, *z3, *z2, col));
	++number_rear;     
      } 
      z1 = z3++;
      if (z3 != path.end()) { 
        z2 = z3++;
      }
    }
    //std::cerr << "... resolved with " << nr << " 2-breaks" << std::endl;
  }
  return number_rear;
}

#endif
