#ifndef STATISTICS_HPP
#define STATISTICS_HPP

template<class mcolor_t>
struct GraphPack<mcolor_t>::Statistics {
  using mcolor_type = mcolor_t;
  
  void calculate(GraphPack<mcolor_t> & graph_pack) { 
    clear();
    count_vertex_statistics(graph_pack);
    count_complete_edges(graph_pack);
    count_connected_components(graph_pack);
    count_rearrangement_statistics(graph_pack);
    count_cycles_statistics(graph_pack);
    count_indel_statistics(graph_pack);
  }
  
  void calculate_history_statistics(GraphPack<mcolor_t> & graph_pack);

private:
  void count_vertex_statistics(GraphPack<mcolor_t> const & graph_pack);
  void count_complete_edges(GraphPack<mcolor_t> const & graph_pack);
  void count_connected_components(GraphPack<mcolor_t> const & graph_pack);
  void count_rearrangement_statistics(GraphPack<mcolor_t> & graph_pack);
  void count_cycles_statistics(GraphPack<mcolor_t> & graph_pack);
  void count_indel_statistics(GraphPack<mcolor_t> & graph_pack);

  void clear() { 
    vertex_statistics.fill(0);
    complete_edges.clear();
    size_component_to_count.clear();
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
  std::array<size_t, 4> vertex_statistics; 
  std::vector<std::pair<vertex_t, vertex_t> > complete_edges; 
  std::map<size_t, size_t> size_component_to_count;

  //vertices
  std::map<mcolor_t, size_t> simple_vertices_count;  	 // simple_vertices_count[min(S,!S)] = # simple vertices incident to S-colored
  std::map<mcolor_t, size_t> simple_vertices_alone_count; // simple_vertices_alone_count[min(S,!S)] = # simple vertices incident to S-colored, with no good neighbors

  //edges
  std::map<std::pair<mcolor_t, mcolor_t>, size_t> multiedges_count; 	// multiedges_count[S] = # multiedges of multicolor S.
	std::map<mcolor_t, size_t> irrer_multiedges_count;	// ME[S] = # irregular multiedges of multicolor S.
 	std::map<mcolor_t, size_t> simple_multiedges_count;	// ME[S] = # simple multiedges of multicolor S.
	
  //cycles
  std::map<mcolor_t, size_t> simple_cycle_count; // cycle of simple vertices
  std::map<mcolor_t, size_t> special_cycle_count; // cycle of simple vertices and oo, of even length

  //indels
  std::vector<std::pair<std::pair<mcolor_t, mcolor_t>, size_t> > complement_indel_stats;

  //history statistics 
  std::map<mcolor_t, size_t> number_twobreaks; 
  std::map<mcolor_t, size_t> number_insertions; 
  std::map<mcolor_t, size_t> number_deletions; 
};

/**
 * Function calculate how many graph have different vertices
 * vertex_statistics[0] - number of vertices
 * vertex_statistics[1] - number of duplication vertices
 * vertex_statistics[2] - number of insertion/deletion vertices
 * vertex_statistics[3] - number of vertices with self loop
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_vertex_statistics(GraphPack<mcolor_t> const & graph_pack) {
  vertex_statistics[0] = graph_pack.graph.size(); 

  for (vertex_t const & x : graph_pack.graph) {
    if (graph_pack.is_duplication_vertex(x)) {
      ++vertex_statistics[1];
    } else if (graph_pack.is_indel_vertex(x)) {
      ++vertex_statistics[2];
    } 

    if (graph_pack.is_have_self_loop(x)) {
      ++vertex_statistics[3];
    }  
  } 
}

/**
 * 
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_complete_edges(GraphPack<mcolor_t> const & graph_pack) { 
  std::unordered_set<vertex_t> processed; 

  for (vertex_t const & x : graph_pack.graph) {
    if (processed.count(x) == 0) {
      mularcs_t const & mularcs = graph_pack.get_all_adjacent_multiedges(x);
      if (mularcs.size() == 1 && mularcs.union_multicolors() == graph_pack.multicolors.get_complete_color()) { 
        complete_edges.push_back(edge_t(x, mularcs.cbegin()->first));
        processed.insert(x); processed.insert(mularcs.cbegin()->first);
      } 
    } 
  }
}

/**
 *
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_connected_components(GraphPack<mcolor_t> const & graph_pack) { 
  utility::equivalence<vertex_t> components = graph_pack.split_on_components(false); 
  using component_t = std::set<vertex_t>;
  std::map<vertex_t, component_t> const & classes = components.get_eclasses<component_t>();
  
  for (auto const &component : classes) {
    ++size_component_to_count[component.second.size()];
  }
}

/*
 *
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_rearrangement_statistics(GraphPack<mcolor_t> & graph_pack) {
  std::unordered_set<vertex_t> processed;

  auto pack_mcolor_lambda = [&] (mcolor_t const & f, mcolor_t const & s) -> std::pair<mcolor_t, mcolor_t> { 
    return std::make_pair(std::max(f, s), std::min(f, s));
  }; 

  for (vertex_t const & x : graph_pack.graph) {
    mularcs_t const & current = graph_pack.get_all_adjacent_multiedges(x); 

    if (graph_pack.is_simple_vertex(x)) {  //we define simple vertices as a regular vertex of multidegree 2. 
      processed.insert(x);
      ++simple_vertices_count[std::min(current.cbegin()->second, current.crbegin()->second)]; 
    }

    for (arc_t const & arc : current) {
    	// count two times, because same underected edge (u, v) and (v, u)
      ++multiedges_count[pack_mcolor_lambda(arc.second, graph_pack.multicolors.get_complement_color(arc.second))]; 

      if ((processed.count(x) != 0) && (processed.count(arc.first) != 0)) {  
				++simple_multiedges_count[arc.second]; //if two vertices have degree = 2 - is simple edges
      } 

      if (arc.first == Infty) { 
        ++multiedges_count[pack_mcolor_lambda(arc.second, graph_pack.multicolors.get_complement_color(arc.second))];
				++multiedges_count[pack_mcolor_lambda(arc.second, graph_pack.multicolors.get_complement_color(arc.second))];
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

/**
 *
 */
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

/**
 *
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::count_indel_statistics(GraphPack<mcolor_t> & graph_pack) { 
  std::map<mcolor_t, size_t> indel_stats;
  std::unordered_set<vertex_t > processed; 

  for (vertex_t const & a1 : graph_pack.graph) {  
    vertex_t const & a2 = graph_pack.graph.get_obverse_vertex(a1);
    if ((processed.count(a1) == 0) && (processed.count(a2) == 0) 
        && graph_pack.is_indel_vertex(a1) && graph_pack.is_indel_vertex(a2))  {
      processed.insert(a1); processed.insert(a2);
    
      mcolor_t const & indel_color = graph_pack.get_all_adjacent_multiedges(a1).union_multicolors(); 
      mcolor_t const & bar_indel_color = graph_pack.multicolors.get_complement_color(indel_color);
      assert(indel_color == graph_pack.get_all_adjacent_multiedges(a2).union_multicolors());
      indel_stats[std::min(bar_indel_color, indel_color)] += 1;
    }
  }

  for (auto const & elem : indel_stats) { 
    complement_indel_stats.push_back(std::make_pair(std::make_pair(elem.first, 
            graph_pack.multicolors.get_complement_color(elem.first)), 
          elem.second));
  }

  using elem_t = std::pair<std::pair<mcolor_t, mcolor_t>, size_t>;
  std::sort(complement_indel_stats.begin(), complement_indel_stats.end(), [](elem_t const & first, elem_t const & second) {
    return (first.second > second.second);
  });
} 

/**
 *
 */
template<class mcolor_t>
void GraphPack<mcolor_t>::Statistics::calculate_history_statistics(GraphPack<mcolor_t> & graph_pack) { 
  for (auto const & twobreak : graph_pack.history) {
    vertex_t const & p = twobreak.get_vertex(0);
    vertex_t const & q = twobreak.get_vertex(1);
    vertex_t const & x = twobreak.get_vertex(2);
    vertex_t const & y = twobreak.get_vertex(3);

    if (p == Infty || q == Infty || x == Infty || y == Infty) {
      ++number_twobreaks[twobreak.get_mcolor()];
    } else { 
      if (p == graph_pack.graph.get_obverse_vertex(x)) { 
        ++number_deletions[twobreak.get_mcolor()];
      } else if (q == graph_pack.graph.get_obverse_vertex(y)) {  
        ++number_deletions[twobreak.get_mcolor()];
      } else if (p == graph_pack.graph.get_obverse_vertex(q)) {  
        ++number_insertions[twobreak.get_mcolor()];
      } else if (x == graph_pack.graph.get_obverse_vertex(y)) { 
        ++number_insertions[twobreak.get_mcolor()];
      } else { 
        ++number_twobreaks[twobreak.get_mcolor()];
      }
    }
  }
}

#endif