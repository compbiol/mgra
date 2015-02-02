#ifndef STATISTICS_HPP
#define STATISTICS_HPP

template<class mcolor_t>
struct GraphPack<mcolor_t>::Statistics {

  void calculate(GraphPack<mcolor_t> & graph_pack) { 
    clear();
    count_vertex_statistics(graph_pack);
    count_indel_statistics(graph_pack);
    count_rearrangement_statistics(graph_pack);
    count_cycles_statistics(graph_pack);
  }
  
private:
  void count_vertex_statistics(GraphPack<mcolor_t> const & graph_pack);
  void count_indel_statistics(GraphPack<mcolor_t> & graph_pack);
  void count_rearrangement_statistics(GraphPack<mcolor_t> const & graph_pack);
  void count_cycles_statistics(GraphPack<mcolor_t> & graph_pack);

  void clear() { 
    vertex_statistics.fill(0);
    simple_vertices_count.clear(); 
    simple_vertices_alone_count.clear(); 
    multiedges_count.clear();
    irrer_multiedges_count.clear();
    simple_multiedges_count.clear();
    simple_cycle_count.clear();
    special_cycle_count.clear();
    complement_indel_stats.clear();
  }
public:
	//summary stat
  std::array<size_t, 3> vertex_statistics; 
  
  //vertices
  std::map<mcolor_t, size_t> simple_vertices_count;  	 // simple_vertices_count[min(S,!S)] = # simple vertices incident to S-colored
  std::map<mcolor_t, size_t> simple_vertices_alone_count; // simple_vertices_alone_count[min(S,!S)] = # simple vertices incident to S-colored, with no good neighbors

  //edges
  std::map<mcolor_t, size_t> multiedges_count; 	// multiedges_count[S] = # multiedges of multicolor S.
	std::map<mcolor_t, size_t> irrer_multiedges_count;	// ME[S] = # irregular multiedges of multicolor S.
 	std::map<mcolor_t, size_t> simple_multiedges_count;	// ME[S] = # simple multiedges of multicolor S.
	
  //cycles
  std::map<mcolor_t, size_t> simple_cycle_count; // cycle of simple vertices
  std::map<mcolor_t, size_t> special_cycle_count; // cycle of simple vertices and oo, of even length

  //indels
  std::map<mcolor_t, size_t> complement_indel_stats;
};

/*
 * Function calculate how many graph have different vertices
 * vertex_statistics[0] - number of duplication vertices
 * vertex_statistics[1] - number of insertion/deletion vertices
 * vertex_statistics[2] - number of vertices with self loop
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_vertex_statistics(GraphPack<mcolor_t> const & graph_pack) {
  for(vertex_t const &x : graph_pack.graph) {
    if (graph_pack.is_duplication_vertex(x)) {
      ++vertex_statistics[0];
    } else if (graph_pack.is_indel_vertex(x)) {
      ++vertex_statistics[1];
    } 

    if (graph_pack.is_have_self_loop(x)) {
      ++vertex_statistics[2];
    }  
  } 
}

template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_indel_statistics(GraphPack<mcolor_t> & graph_pack) { 
  std::unordered_set<vertex_t > processed; 

  for (vertex_t const & a1 : graph_pack.graph) {  
    vertex_t const & a2 = graph_pack.graph.get_obverse_vertex(a1);
    if ((processed.count(a1) == 0) && (processed.count(a2) == 0) 
    		&& graph_pack.is_indel_vertex(a1) && graph_pack.is_indel_vertex(a2))  {
      processed.insert(a1); processed.insert(a2);
    
      mcolor_t const & indel_color = graph_pack.get_all_adjacent_multiedges(a1).union_multicolors(); 
      mcolor_t const & bar_indel_color = graph_pack.multicolors.get_complement_color(indel_color);
      assert(indel_color == graph_pack.get_all_adjacent_multiedges(a2).union_multicolors());
      complement_indel_stats[std::min(bar_indel_color, indel_color)] += 1;
    }
  }

} 

template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_rearrangement_statistics(GraphPack<mcolor_t> const & graph_pack) {
  std::unordered_set<vertex_t> processed;

  for (vertex_t const & x : graph_pack.graph) {
    mularcs_t const & current = graph_pack.get_all_adjacent_multiedges(x); 

    if (graph_pack.is_simple_vertex(x)) {  //we define simple vertices as a regular vertex of multidegree 2. 
      processed.insert(x);
      ++simple_vertices_count[std::min(current.cbegin()->second, current.crbegin()->second)]; 
    }

    for (arc_t const & arc : current) {
    	// count two times, because same underected edge (u, v) and (v, u)
      ++multiedges_count[arc.second]; 

      if ((processed.count(x) != 0) && (processed.count(arc.first) != 0)) {  
				++simple_multiedges_count[arc.second]; //if two vertices have degree = 2 - is simple edges
      } 

      if (arc.first == Infty) { 
				++multiedges_count[arc.second];
				++irrer_multiedges_count[arc.second];			
      } 
    }
  }

  // count lonely vertices (short paths) 
  for (vertex_t const & v : processed) {
    mularcs_t const & current = graph_pack.get_all_adjacent_multiedges(v);
    if (processed.find(current.cbegin()->first) == processed.end() && processed.find(current.crbegin()->first) == processed.end()) {
      ++simple_vertices_alone_count[std::min(current.cbegin()->second, current.crbegin()->second)]; //no good neighbors
    }
  } 
} 

template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_cycles_statistics(GraphPack<mcolor_t> & graph_pack) { 
  std::unordered_set<vertex_t> processed;

  for(vertex_t const & x : graph_pack.graph) { 
    if (processed.count(x) != 0) continue;  

    mularcs_t const & mularcs_x = graph_pack.get_all_adjacent_multiedges(x); 
    if (!(graph_pack.is_simple_vertex(x) 
      && graph_pack.multicolors.get_complement_color(mularcs_x.cbegin()->second) == mularcs_x.crbegin()->second)) { 
      continue;
    } 

    vertex_t current = x;
    vertex_t prev = "";
    mcolor_t special_Q; 

    do {
      processed.insert(current);

      if (!graph_pack.is_simple_vertex(current)) break;

      mularcs_t const & mularcs_y = graph_pack.get_all_adjacent_multiedges(current);

      if (prev == mularcs_y.cbegin()->first) {
        prev = current;
        current = mularcs_y.crbegin()->first;
      } else {
        prev = current;
        current = mularcs_y.cbegin()->first;
      }

      while (current == Infty) {
      	if (special_Q.empty()) {
      	  special_Q = graph_pack.get_all_multicolor_edge(prev, current); 
      	  prev = x;
      	  current = mularcs_x.cbegin()->first; 
      	} else {
      	  if (special_Q != graph_pack.get_all_multicolor_edge(prev, current)) {
      	    ++special_cycle_count[std::min(special_Q, graph_pack.get_all_multicolor_edge(prev, current))]; 	  
      	  }
      	  break;
      	}
      }
    } while ((current != Infty) && (processed.count(current) == 0));
	
    if (current == x) { //find cycle. 
      ++simple_cycle_count[std::min(mularcs_x.cbegin()->second, mularcs_x.crbegin()->second)];
    }
  }
} 

#endif