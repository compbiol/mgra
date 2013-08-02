#ifndef STAGE1_H_
#define STAGE1_H_

template<class graph_t>
bool Algorithm<graph_t>::stage1() {
  bool symplified = false; 
  size_t num_rear = 0; // number of rearrangements 

  do {
    num_rear = 0; 
    for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
      if (!graph.is_simple_vertex(*is)) { 
	continue; 
      } 

      Mularcs<Mcolor> current = graph.get_adjacent_multiedges(*is);

      path_t path({*is});

      std::unordered_set<vertex_t> processed({*is, Infty}); // we count oo as already processed

      for(auto im = current.cbegin(); im != current.cend(); ++im) {	
	bool is_next = (im == current.cbegin()); 

	vertex_t current = find_simple_path(path, processed, *is, im->first, is_next);

	if (current == *is) { 
	  break; // got a cycle from x to x, cannot extend it 
	}  		    
      }

      if ((*path.begin() == *path.rbegin()) && (*path.begin() != Infty) 
	&& (graph.is_duplication_vertex(*path.begin()) || graph.is_indel_vertex(*path.begin())) 
	&& (path.size() % 2 == 0)) {
	continue; 
      } else {
     	num_rear += process_simple_path(path);
      }
    } 

    if (num_rear != 0) { 
      symplified = true;
    } 
  } while (num_rear > 0); 
   
  return symplified;
} 

template<class graph_t>
vertex_t Algorithm<graph_t>::find_simple_path(path_t& path, std::unordered_set<vertex_t>& processed, const vertex_t& prev, const vertex_t& cur, bool is_next) { 
  std::string previous  = prev;
  std::string current = cur;

  while (true) {
    //FIXME: is_duplication_vertice work is a long while. And uses iff prevent duplication vertex. x -> ... -> y -> z -> t , 
    //if z - end path and edge colors z->t, y->z  complimentary, but t - is duplication vertex and 2-break down all colors.       

    if (is_next) { 
      path.push_front(current);
    } else { 
      path.push_back(current);
    } 

    if (processed.find(current) != processed.end()) { 
      break;
    }
 
    processed.insert(current);

    if (current != Infty && (graph.is_duplication_vertex(current) || graph.is_indel_vertex(current))) { 
	break;
    } 

    Mularcs<Mcolor> new_edges = graph.get_adjacent_multiedges(current);
    Mcolor previous_color = new_edges.find(previous)->second;
    new_edges.erase(previous);
    
    if (new_edges.size() != 1 
    	|| !graph.is_T_consistent_color(new_edges.cbegin()->second) 
	|| !graph.is_T_consistent_color(previous_color) 
	|| graph.get_complement_color(previous_color) != new_edges.cbegin()->second) 
    { 
      break;
    } 

    previous = current;
    current = new_edges.cbegin()->first;
  }

  return current;
} 

template<class graph_t>
size_t Algorithm<graph_t>::process_simple_path(path_t& path) {
  size_t nr = 0;

  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
    //std::cerr << std::endl << "Processing a path of length " << path.size() - 1 << std::endl;
    //std::cerr << "path:\t" << *path.begin();
    //for(auto ip = ++path.begin(); ip != path.end(); ++ip) {
      //std::cerr << " -- " << *ip;
    //}
    //std::cerr << std::endl;

    if (path.size() % 2 && (*path.begin() != *path.rbegin())) {
      //std::cerr << "... ";
      if (!graph.is_vec_T_color(graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin())) ) { //FIXME
	path.erase(path.begin());
	//std::cerr << "left";
      } else {
	path.erase(--path.end());
	//std::cerr << "right";
      }
      //std::cerr << " end removed" << std::endl;
    }

    /*if (*path.begin() == Infty && *path.rbegin() == Infty ) {
      outlog << "... affecting two chromosome ends" << std::endl;
    } else if( *path.begin() == Infty || *path.rbegin() == Infty ) {
      outlog << "... affecting a chromosome end" << std::endl;
    }*/

    if (*path.begin() == *path.rbegin()) {
      if (path.size() % 2 == 0) {
	if (*path.begin() != Infty) {
	  //std::cerr << "ERROR: Semi-cycle w/o infinity! " << *path.begin() << std::endl;
	  exit(1);
	}
	if (graph.is_vec_T_color(graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin()))) { //FIXME
	  //std::cerr << "... semi-cycle, fusion applied" << std::endl;
	  const vertex_t& x0 = *(++path.begin());
	  const vertex_t& y0 = *(++path.rbegin());
    
	  TwoBreak<Mcolor> t(Infty, x0, Infty, y0, graph.get_adjacent_multiedges(x0).get_multicolor(Infty)); 

	  graph.apply_two_break(t);
	  path.erase(--path.end());
	  *path.begin() = y0;
	  ++nr;
	} else {
	  //std::cerr << "... semi-cycle, fission applied" << std::endl;
	  const std::string y0 = *(++path.rbegin());
	  const std::string y1 = *(++++path.rbegin());

	  graph.apply_two_break(TwoBreak<Mcolor>(y0, y1, Infty, Infty, graph.get_adjacent_multiedges(y0).get_multicolor(y1)));
	  ++nr;
	  path.erase(--path.end());
	  *path.rbegin() = Infty;
	}
	if (path.size() < 4) return nr;
      } else { 
	;//std::cerr << "... cycle" << std::endl;
      } 
    }

    Mcolor Q;

    while (true) {
      // multicolor of (z1,z2). N.B.: x2 is NOT oo
      //std::cerr << "... multicolors of first and second multiedges: ";
    
      Q = graph.get_adjacent_multiedges(*(++path.begin())).get_multicolor(*path.begin());
    
      if (graph.is_vec_T_color(Q)) { 
	break;
      } 

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
	    return nr;
	  }
	}
      }
    }

    Mcolor Qrep = graph.get_min_complement_color(Q);

    // x1 -- x2 -- x3 -- ... -- x2k
    // results in:
    // x2 == x3   x4 == x5 ... x(2k-2) == x(2k-1)   x1 -- x2k

    path_t::const_iterator z3 = path.begin();
    path_t::const_iterator z0 = z3++;
    path_t::const_iterator z1 = z3++;
    path_t::const_iterator z2 = z3++;

    while(z3 != path.end()) {
	graph.apply_two_break(TwoBreak<Mcolor>(*z0, *z1, *z3, *z2, Q));
	++nr;
      z1 = z3++;
      if (z3 == path.end()) { 
	break;
      } 
      z2 = z3++;
    }

    //std::cerr << "... resolved with " << nr << " 2-breaks" << std::endl;
    return nr;
  }

  return 0;
}

#endif
