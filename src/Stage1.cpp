#include "Stage1.h"

size_t Stage1::process_simple_path(/*path_t& path, MBGraph& graph*/) {
  size_t nr = 0;

  if (path.size() >= 4 || (path.size() == 3 && *path.begin() == *path.rbegin())) {
    outlog << std::endl << "Processing a path of length " << path.size() - 1 << std::endl;
    outlog << "path:\t" << *path.begin();
    for(auto ip = ++path.begin(); ip != path.end(); ++ip) {
      outlog << " -- " << *ip;
    }
    outlog << std::endl;

    if (path.size() % 2 && (*path.begin() != *path.rbegin())) {
      outlog << "... ";
      if (!member(graph.DiColor, graph.get_adjacent_multiedges(*(++path.begin()))[*path.begin()] ) ) {
	path.erase(path.begin());
	outlog << "left";
      } else {
	path.erase(--path.end());
	outlog << "right";
      }
      outlog << " end removed" << std::endl;
    }

    if (*path.begin() == Infty && *path.rbegin() == Infty ) {
      outlog << "... affecting two chromosome ends" << std::endl;
    } else if( *path.begin() == Infty || *path.rbegin() == Infty ) {
      outlog << "... affecting a chromosome end" << std::endl;
    }

    if (*path.begin() == *path.rbegin()) {
      if (path.size() % 2 == 0) {
	if (*path.begin() != Infty) {
	  outlog << "ERROR: Semi-cycle w/o infinity!" << std::endl;
	  exit(1);
	}
	if (member(graph.DiColor, graph.get_adjacent_multiedges(*(++path.begin()))[*path.begin()]) ) {
	  outlog << "... semi-cycle, fusion applied" << std::endl;
	  const std::string& x0 = *(++path.begin());
	  const std::string& y0 = *(++path.rbegin());
    
	  TwoBreak t(Infty, x0, Infty, y0, graph.get_adjacent_multiedges(x0)[Infty]);
	  if(t.apply(graph, true)) {
	    path.erase(--path.end());
	    *path.begin() = y0;
	    ++nr;
	  }
	} else {
	  outlog << "... semi-cycle, fission applied" << endl;
	  const std::string y0 = *(++path.rbegin());
	  const std::string y1 = *(++++path.rbegin());

	  if (TwoBreak(y0, y1, Infty, Infty, graph.get_adjacent_multiedges(y0)[y1]).apply(graph, true)) {
	    ++nr;
	    path.erase(--path.end());
	    *path.rbegin() = Infty;
	  }
	}
	if (path.size() < 4) return nr;
      } else { 
	outlog << "... cycle" << std::endl;
      } 
    }

    Mcolor Q;

    while (true) {
      // multicolor of (z1,z2). N.B.: x2 is NOT oo
      outlog << "... multicolors of first and second multiedges: ";
    
      Q = graph.get_adjacent_multiedges(*(++path.begin()))[*path.begin()];
    
      outlog << genome_match::mcolor_to_name(Q) << " + " << genome_match::mcolor_to_name(graph.CColor(Q)) << endl;

      if (member(graph.DiColor, Q)) { 
	break;
      } 

      if (*path.begin() == *path.rbegin()) {
	outlog << "... rotating" << std::endl;
	path.push_back(*path.begin());
	path.erase(path.begin());
      } else {
	if (*path.begin() == Infty && *path.rbegin() != Infty) {
	  outlog << "... flipping" << std::endl;
	  for(auto ip = ++path.begin();ip != path.end();) {
	    path.push_front(*ip);
	    path.erase(ip++);
	  }
	}
	if (*path.rbegin() == Infty) {
	  outlog << "... extending beyond oo" << std::endl;
	  path.push_back(Infty);
	  path.erase(path.begin());
	} else {
	  outlog << "... truncating ??" << std::endl;
	  path.erase(path.begin());
	  path.erase(--path.end());
	  if (path.size() < 4) return nr;
	}
      }
    }

    Mcolor Qrep = graph.CColorRep(Q);

    // x1 -- x2 -- x3 -- ... -- x2k
    // results in:
    // x2 == x3   x4 == x5 ... x(2k-2) == x(2k-1)   x1 -- x2k

    path_t::const_iterator z3 = path.begin();
    path_t::const_iterator z0 = z3++;
    path_t::const_iterator z1 = z3++;
    path_t::const_iterator z2 = z3++;

    while(z3 != path.end()) {
      if( TwoBreak(*z0, *z1, *z3, *z2, Q).apply(graph, true) ) {
	++nr;
      } else {
	z0 = z2;
      }
      z1 = z3++;
      if (z3 == path.end()) { 
	break;
      } 
      z2 = z3++;
    }

    outlog << "... resolved with " << nr << " 2-breaks" << std::endl;
    return nr;
  }

  return 0;
}

vertex_t Stage1::find_simple_path(/*path_t& path, MBGraph& graph, std::unordered_set<vertex_t>& processed,*/ const vertex_t& prev, const vertex_t& cur, bool is_next) { 
  std::string previous  = prev;
  std::string current = cur;

  while (true) {
    if (is_next) { 
      path.push_front(current);
    } else { 
      path.push_back(current);
    } 

    if (processed.find(current) != processed.end()) { 
      break;
    } 
    processed.insert(current);

    mularcs_t new_edges = graph.get_adjacent_multiedges(current);
    new_edges.erase(previous);

    if (!(new_edges.size() == 1 && new_edges.begin()->second.is_good_multiedge()) || !graph.is_T_consistent_color(new_edges.begin()->second)) { 
      break;
    } 

    previous = current;
    current = new_edges.begin()->first;
  }

  return current;
} 

bool Stage1::stage1(/*MBGraph& graph*/) {
  bool symplified = false; 
  size_t nr = 0; // number of rearrangements 

  do {
    nr = 0; 
    for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
      //multimularcs_t current = graph.get_adjacent_multiedges_v2(*is);
      mularcs_t current = graph.get_adjacent_multiedges(*is);

      if (!graph.is_simple_vertice(current)) { 
	continue; 
      } 

      path.clear(); 
      path.push_back(*is);

      processed.clear(); 
      processed.insert(*is); 
      processed.insert(Infty); // we count oo as already processed

      for(auto im = current.begin(); im != current.end(); ++im) {
	if (!graph.is_T_consistent_color(im->second)) { 
	  continue; // not T-consistent
	} 
	
	bool is_next = (im == current.begin()); 
	std::string current = find_simple_path(/*path, graph, processed,*/ *is, im->first, is_next);

	if (current == *is) { 
	  break; // got a cycle from x to x, cannot extend it 
	}  		    
      }
      nr += process_simple_path(/*path, graph*/);
    } 

    if (nr != 0) { 
      symplified = true;
    } 
  } while (nr > 0); 
   
  return symplified;
} 
