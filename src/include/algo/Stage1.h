#ifndef STAGE1_H_
#define STAGE1_H_

#include <list>
#include <set>
#include <string>

#include "mpbgraph.h"
#include "2break.h"

typedef std::list<vertex_t> path_t;

//add graph - shared_ptr

//Stage 1: loop over vertices  
template<class graph_t>
struct Stage1 {  
	Stage1(graph_t& gr)
	: graph(gr) { 
	} 

	bool stage1(); 
	
	graph_t get_graph() {
		return graph;
 	}
private: 
	size_t process_simple_path();	
	vertex_t find_simple_path(const vertex_t& prev, const vertex_t& cur, bool is_next); 
private: 
	path_t path;
	std::unordered_set<vertex_t> processed;
	graph_t graph;
};

template<class graph_t>
bool Stage1<graph_t>::stage1() {
  bool symplified = false; 
  size_t num_rear = 0; // number of rearrangements 

  do {
    num_rear = 0; 
    for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
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
	std::string current = find_simple_path(*is, im->first, is_next);

	if (current == *is) { 
	  break; // got a cycle from x to x, cannot extend it 
	}  		    
      }
      num_rear += process_simple_path();
    } 

    if (num_rear != 0) { 
      symplified = true;
    } 
  } while (num_rear > 0); 
   
  return symplified;
} 

template<class graph_t>
vertex_t Stage1<graph_t>::find_simple_path(const vertex_t& prev, const vertex_t& cur, bool is_next) { 
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

template<class graph_t>
size_t Stage1<graph_t>::process_simple_path() {
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

#endif
