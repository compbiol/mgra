#ifndef STAGE1_H_
#define STAGE1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage1() {
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  do {
    number_rear = 0; 
    for (auto const & v: *graph) {  
      if (graph->is_simple_vertex(v)) { 
	path_t path({v});
	std::unordered_set<vertex_t> processed({v, Infty}); // we count oo as already processed

        auto const  find_simple_path_lambda = [&] (vertex_t const & prev, vertex_t const & cur, bool is_next) -> vertex_t { 
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

            if (processed.find(current) == processed.end() && !graph->is_duplication_vertex(current)) {     
              processed.insert(current);
              mularcs_t&& new_edges = graph->get_adjacent_multiedges(current);
              auto const & previous_color = new_edges.get_multicolor(previous); 
              new_edges.erase(previous);
    
              if (new_edges.size() == 1 && graph->get_complement_color(previous_color) == new_edges.cbegin()->second) {
	        mularcs_t const & edges = graph->get_adjacent_multiedges_with_info(current, false);
                auto const count_Lambda = [&] (const vertex_t& v) -> bool {
                  bool flag = true; 
                  auto const & colors = edges.equal_range(v);
                  for (auto arc = colors.first; arc != colors.second && flag; ++arc) {
                    flag = graph->is_vec_T_consistent_color(arc->second);
                  }     
                  return flag;
                };

                bool nedge = count_Lambda(new_edges.cbegin()->first); 
                bool pedge = count_Lambda(previous);
                if (nedge || pedge) { 
                  previous = current;
	          current = new_edges.cbegin()->first; 
	          stop = true;
                }
              } 
            }
          } 
          return current;
        };

        mularcs_t const & current = graph->get_adjacent_multiedges(v);
        //std::cerr << std::endl << "start " << v << std::endl; 
        //std::cerr << "go to " << current.cbegin()->first << genome_match::mcolor_to_name(current.cbegin()->second) << std::endl; 
        vertex_t const & last = find_simple_path_lambda(v, current.cbegin()->first, true);
        if (last != v) { 
          //std::cerr << "go to " << current.crbegin()->first << genome_match::mcolor_to_name(current.crbegin()->second) << std::endl; 
          find_simple_path_lambda(v, current.crbegin()->first, false);
        }  		    

        number_rear += process_simple_path(path);
      }
    } 

    if (number_rear != 0) { 
      isChanged = true;
    } 
  } while (number_rear > 0); 
   
  return isChanged;
} 

template<class graph_t>
size_t Algorithm<graph_t>::process_simple_path(path_t& path) {
  size_t number_rear = 0;

  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
#ifdef LOG_ENABLED
    std::cerr << std::endl << "Processing a path of length " << path.size() - 1 << std::endl;
    std::cerr << "path:\t" << *path.begin();
    for(auto ip = ++path.begin(); ip != path.end(); ++ip) {
      std::cerr << " -- " << *ip;
    }
    std::cerr << std::endl;
#endif
    auto const & edges = graph->get_adjacent_multiedges_with_info(*(++path.begin()), false);
    auto const count_Lambda = [&] (vertex_t const & v) -> bool {
      bool flag = true; 
      auto const & colors = edges.equal_range(v);
      for (auto arc = colors.first; arc != colors.second && flag; ++arc) {
        flag = graph->is_vec_T_consistent_color(arc->second);
      }     
      return flag;
    };

    mcolor_t process_color; 
    bool pedge = count_Lambda(*path.begin());
    bool nedge = count_Lambda(*(++++path.begin())); 
    std::set<mcolor_t> process_colors;

    if (pedge) { 
      process_color = graph->get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
      auto const & mul = graph->get_adjacent_multiedges_with_info(*(++path.begin()), false);
      auto const & pair_colors = mul.equal_range(*path.begin());
      for (auto col = pair_colors.first; col != pair_colors.second; ++col) {
 	process_colors.insert(col->second);     
      } 
    } else if (!pedge && nedge) { 
      process_color = graph->get_adjacent_multiedges(*(++path.begin())).get_multicolor(*(++++path.begin()));
      auto const & mul = graph->get_adjacent_multiedges_with_info(*(++path.begin()), false);
      auto const & pair_colors = mul.equal_range(*(++++path.begin())); 
      for (auto col = pair_colors.first; col != pair_colors.second; ++col) {
       	process_colors.insert(col->second);
      } 
    } 
    
    if ((path.size() % 2 != 0) && (*path.begin() != *path.rbegin())) {
      if (process_color != graph->get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin())) { 
	path.erase(path.begin());
      } else {
	path.erase(--path.end());
      }
    }

    if (*path.begin() == *path.rbegin()) {
      if (path.size() % 2 == 0) {
	if (process_color == graph->get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin())) { 
#ifdef LOG_ENABLED
	  std::cerr << "... semi-cycle, fusion applied" << std::endl;
#endif
	  vertex_t const & self_v = *(path.begin());
	  vertex_t const & x0 = *(++path.begin());
	  vertex_t const & y0 = *(++path.rbegin());

	  mularcs_t const & mul = graph->get_adjacent_multiedges_with_info(x0, false);
	  auto const & colors = mul.equal_range(self_v);
	  for (auto it = colors.first; it != colors.second; ++it) { 
	    graph->apply_two_break(twobreak_t(self_v, x0, self_v, y0, it->second));
	    ++number_rear;
	  } 
	
	  path.erase(--path.end());
	  *path.begin() = y0;
	} else {
#ifdef LOG_ENABLED
	  std::cerr << "... semi-cycle, fission applied" << std::endl;
#endif
	  vertex_t const & self_v = *(path.begin());
	  vertex_t const & y0 = *(++path.rbegin());
	  vertex_t const & y1 = *(++++path.rbegin());

	  mularcs_t const & mul = graph->get_adjacent_multiedges_with_info(y0, false);
	  auto const & pair = mul.equal_range(y1);
	  for (auto it = pair.first; it != pair.second; ++it) { 
 	    graph->apply_two_break(twobreak_t(y0, y1, self_v, self_v, it->second));
	    ++number_rear;
	  }
        
	  path.erase(--path.end());
	  *path.rbegin() = self_v;
	}

	if (path.size() < 4) { 
	  return number_rear;
	}
      } 
#ifdef LOG_ENABLED
      else { 
	std::cerr << "... cycle" << std::endl;
      }
#endif 
    }
  
    auto color = graph->get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
    while (process_color != color) {
      // multicolor of (z1,z2). N.B.: x2 is NOT oo
#ifdef LOG_ENABLED
      std::cerr << "... multicolors of first and second multiedges: ";
#endif
      if (*path.begin() == *path.rbegin()) {
#ifdef LOG_ENABLED
	std::cerr << "... rotating" << std::endl;
#endif
	path.push_back(*path.begin());
	path.erase(path.begin());
      } else {
	if (*path.begin() == Infty && *path.rbegin() != Infty) {
#ifdef LOG_ENABLED
	  std::cerr << "... flipping" << std::endl;
#endif
	  for(auto ip = ++path.begin(); ip != path.end();) {
	    path.push_front(*ip);
	    path.erase(ip++);
	  }
	}
	if (*path.rbegin() == Infty) {
#ifdef LOG_ENABLED
	  std::cerr << "... extending beyond oo" << std::endl;
#endif
	  path.push_back(Infty);
	  path.erase(path.begin());
	} else {
#ifdef LOG_ENABLED
	  std::cerr << "... truncating ??" << std::endl;
#endif
	  path.erase(path.begin());
	  path.erase(--path.end());
	  if (path.size() < 4) { 
	    return number_rear;
	  }
	}
      }

      color = graph->get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
    }

    // x1 -- x2 -- x3 -- ... -- x2k
    // results in:
    // x2 == x3   x4 == x5 ... x(2k-2) == x(2k-1)   x1 -- x2k

    //std::cerr << genome_match::mcolor_to_name(process_color) << std::endl;
     
    auto z3 = path.begin();
    auto z0 = z3++;
    auto z1 = z3++;
    auto z2 = z3++;
    
    while (z3 != path.end()) {
      for (auto const & col : process_colors) {
	graph->apply_two_break(twobreak_t(*z0, *z1, *z3, *z2, col));
	++number_rear;     
      } 
      z1 = z3++;
      if (z3 != path.end()) { 
        z2 = z3++;
      }
    }

#ifdef LOG_ENABLED
    std::cerr << "... resolved with " << number_rear << " 2-breaks" << std::endl;
#endif
  }
  return number_rear;
}

#endif
