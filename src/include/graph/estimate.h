#ifndef ESTIMATE_H_
#define ESTIMATE_H_

template<class graph_t>
struct Statistics { 
  typedef structure::Mcolor mcolor_t; 
  typedef structure::Mularcs<mcolor_t> mularcs_t; 

  Statistics(const std::shared_ptr<graph_t>& gr)
  : graph(gr) 
  {
    count_vertex_statistics();   
    count_compl_multiedges();
    count_cycles();
    count_indel_statistics();
  }	      	
 
  std::array<size_t, 4> get_vertex_statistics() const {
    return vertex_statistics;
  } 

  std::vector<std::string> get_compl_stat();   

  std::multimap<size_t, std::tuple<mcolor_t, mcolor_t, size_t, size_t> > get_indel_stat() const { 
    return indel_stats; 
  }

private:
  void count_vertex_statistics();
  void count_compl_multiedges(); 
  void count_indel_statistics();
  void count_cycles();

private: 
  std::shared_ptr<graph_t> graph;

  std::array<size_t, 4> vertex_statistics; 
  
  //vertices
  std::unordered_map<size_t, size_t > multidegree_count; // multidegree_count[n] = # vertices of multidegree n. 
  std::map<mcolor_t, size_t> simple_vertices_count;  	 // simple_vertices_count[min(S,!S)] = # simple vertices incident to S-colored
  std::map<mcolor_t, size_t> simple_vertices_alone_count;  // simple_vertices_alone_count[min(S,!S)] = # simple vertices incident to S-colored, with no good neighbors

  //edges
  std::map<mcolor_t, size_t> not_compl_multiedges_count;	// not complement multiedges[S] = # multiedges of not complement multicolor S,
  std::map<mcolor_t, size_t> compl_multiedges_count; 	// multiedges_count[S] = # multiedges of multicolor S.
  std::map<mcolor_t, size_t> good_multiedges_count; 	// good_multiedges_count[S] = # good multiedges of multicolor S. 
  std::map<mcolor_t, size_t> good_irrer_multiedges_count;	// ME[S] = # good irregular multiedges of multicolor S.
  std::map<mcolor_t, size_t> simple_multiedges_count;	// ME[S] = # simple multiedges of multicolor S.
	
  //cycles
  std::map<mcolor_t, size_t> simple_cycle_count; 		// cycle of simple vertices
  std::map<mcolor_t, size_t> special_cycle_count; 	// cycle of simple vertices and oo, of even length

  std::multimap<size_t, std::tuple<mcolor_t, mcolor_t, size_t, size_t> > indel_stats;
};

template<class graph_t>
void Statistics<graph_t>::count_vertex_statistics() {
  vertex_statistics.fill(0);

  for(const auto &x : *graph) {
    if (graph->is_duplication_vertex(x)) {
	++vertex_statistics[0];
    } else if (graph->is_indel_vertex(x)) {
	++vertex_statistics[1];
    } 

    if (graph->is_have_self_loop(x)) {
	++vertex_statistics[2];
    } 

    const mularcs_t& current = graph->get_adjacent_multiedges(x); //current is list with adjacent multiedge&
    for (auto it = current.cbegin(); it != current.cend(); ++it) {
	if (!it->second.is_one_to_one_match()) {
	  ++vertex_statistics[3];
	}
    } 
  } 
} 

template<class graph_t>
void Statistics<graph_t>::count_compl_multiedges() {
  std::unordered_set<std::string> processed;

  for(const auto &x : *graph) {
    const mularcs_t& current = graph->get_adjacent_multiedges(x); //current is list with adjacent multiedges

    ++multidegree_count[current.size()]; //current.size - is degree vertex *it

    if (graph->is_simple_vertex(x)) {  //we define simple vertices as a regular vertex of multidegree 2. 
      processed.insert(x);
      ++simple_vertices_count[std::min(current.cbegin()->second, current.crbegin()->second)]; //simple vertices because degree 2.
    }

    for(auto im = current.cbegin(); im != current.cend(); ++im) {
      ++compl_multiedges_count[im->second];   // count two times, because same underected edge (u, v) and (v, u)
			
      if (graph->is_simple_vertex(x)) {
	++good_multiedges_count[im->second]; //good if one vertices have degree 2
	
	if (im->first == Infty) { 
	  ++good_multiedges_count[im->second];
	} 
						
	if (processed.find(im->first) != processed.end()) {  
	  ++simple_multiedges_count[im->second]; //if two vertices have degree = 2 - is simple edges
	} 
      } 

      if (im->first == Infty) { 
	++compl_multiedges_count[im->second];
	++good_irrer_multiedges_count[im->second];			
      } 
    }
  }

  // count lonely vertices (short paths) 
  for(const auto& v : processed) {
    const mularcs_t& current = graph->get_adjacent_multiedges(v);
    if (processed.find(current.cbegin()->first) == processed.end() && processed.find(current.crbegin()->first) == processed.end()) {
      ++simple_vertices_alone_count[std::min(current.cbegin()->second, current.crbegin()->second)]; //no good neighbors
    }
  } 
} 

template<class graph_t>
void Statistics<graph_t>::count_cycles() { 
  std::unordered_set<vertex_t> processed;

  for(const auto &x : *graph) { 
    if (processed.count(x) != 0) { 
      continue; 
    } 

    const mularcs_t& mularcs_x = graph->get_adjacent_multiedges(x); 

    if (!(graph->is_simple_vertex(x) && graph->get_complement_color(mularcs_x.cbegin()->second) == mularcs_x.crbegin()->second)) { 
      continue;
    } 

    vertex_t current = x;
    vertex_t prev = "";
    mcolor_t special_Q; 

    do {
      processed.insert(current);

      if (!graph->is_simple_vertex(current)) {
	break;
      }

      const mularcs_t& mularcs_y = graph->get_adjacent_multiedges(current);

      if (prev == mularcs_y.cbegin()->first) {
	prev = current;
	current = mularcs_y.crbegin()->first;
      } else {
	prev = current;
	current = mularcs_y.cbegin()->first;
      }

      while (current == Infty) {
	if (special_Q.empty()) {
	  special_Q = mularcs_y.get_multicolor(current);
	  prev = x;
	  current = mularcs_x.cbegin()->first; 
	} else {
	  if (special_Q != mularcs_y.get_multicolor(current)) {
	    ++special_cycle_count[std::min(special_Q, mularcs_y.get_multicolor(current))]; 	  
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

template<class graph_t>
void Statistics<graph_t>::count_indel_statistics() { 
  std::map<std::pair<mcolor_t, mcolor_t>, size_t> temp;
  std::unordered_set<vertex_t > processed; 

  for (const auto &a1 : *graph) {  
    const vertex_t& a2 = graph->get_obverse_vertex(a1);
    const mularcs_t mularcs = graph->get_adjacent_multiedges(a1);

    if (graph->is_indel_vertex(a1) && (processed.count(a1) == 0) && graph->is_indel_vertex(a2))  {
      //std::cerr << "Start worked with " << a1 << " " << a2;
      processed.insert({a1, a2});

      const auto& indel_color = mularcs.union_multicolors(); 
      const auto& bar_indel_color = graph->get_complement_color(indel_color);
      assert(indel_color == graph->get_adjacent_multiedges(a2).union_multicolors());

      if (temp.count(std::make_pair(bar_indel_color, indel_color)) == 0) {
        temp.insert(std::make_pair(std::make_pair(bar_indel_color, indel_color), 1));
      } else { 	
      	temp.find(std::make_pair(bar_indel_color, indel_color))->second += 1;
      }
    }
  }
	
  for (const auto &stat : temp) {
    std::tuple<mcolor_t, mcolor_t, size_t, size_t> infor(stat.first.first, stat.first.second, graph->split_color(stat.first.first, false).size(), graph->split_color(stat.first.second, false).size());
    indel_stats.insert(std::make_pair(stat.second, infor));
  } 
} 

template<class graph_t>
std::vector<std::string> Statistics<graph_t>::get_compl_stat() { 
  auto calc_value = [] (const std::map<mcolor_t, size_t>& where, const mcolor_t& what) -> size_t { 
    if (where.find(what) != where.end()) { 
      return where.find(what)->second;
    } 
    return 0; 
  }; 

  std::multimap<size_t, std::string> answer;

  for(auto im = compl_multiedges_count.cbegin(); im != compl_multiedges_count.cend(); ++im) {
    const auto& current = graph->get_complement_color(im->first);  // complementary multicolor.

    if (im->first < current) {
      continue;
    } 
	 
    size_t m1 = calc_value(compl_multiedges_count, current) / 2;
    size_t m2 = (im->second) / 2;	
    int paths = calc_value(simple_vertices_count, current) - (calc_value(simple_multiedges_count, current) + calc_value(simple_multiedges_count, im->first)) - calc_value(simple_vertices_alone_count, current) - calc_value(special_cycle_count, current);
    int cycles = calc_value(simple_cycle_count, current) + calc_value(special_cycle_count, current);
 	
    std::ostringstream os;
 
    os << "{";		    
    if (graph->is_T_consistent_color(im->first)) { 
      os << "\\bf ";
    } 

    const auto& first = graph->get_min_complement_color(current); 
    const auto& second = graph->get_complement_color(first);

    os <<  genome_match::mcolor_to_name(first) << " + "  <<  genome_match::mcolor_to_name(second) << "} & " 
      // multiedges
       << m1 << " + " << m2 << " = " << m1 + m2 << " & " 
      // simple vertices
       << calc_value(simple_vertices_count, current) << " & "  
      // simple multiedges
       << calc_value(simple_multiedges_count, current) << " + " << calc_value(simple_multiedges_count, im->first) << " = " << calc_value(simple_multiedges_count, current) + calc_value(simple_multiedges_count, im->first) << " & "
      // simple paths + cycles
       << paths << " + " << cycles << " = " << paths + cycles << " & "
      // irregular multiedges
       << calc_value(good_irrer_multiedges_count, current) << " + " << calc_value(good_irrer_multiedges_count, im->first) << " = " << calc_value(good_irrer_multiedges_count, current) + calc_value(good_irrer_multiedges_count, im->first);
	
    answer.insert(std::make_pair(m1 + m2, os.str()));
  }
	
  std::vector<std::string> output;
  for(auto it = answer.rbegin(); it != answer.rend(); ++it) {
    output.push_back(it->second); 
  }  	
  return output;
} 

#endif
