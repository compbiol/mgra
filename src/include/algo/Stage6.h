#ifndef STAGE6_H_ 
#define STAGE6_H_ 

template<class graph_t>
size_t Algorithm<graph_t>::calculate_cost(const vertex_t& y, const Mularcs<Mcolor>& mularcs_x, const Mularcs<Mcolor>& mularcs_y) { 
  typedef std::pair<std::pair<vertex_t, Mcolor>, size_t> colacr_t;
  equivalence<colacr_t> equiv; 

  for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
    if (arc_x->first == y) { 
      continue;
    } 
    for (auto arc_y =  mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
      Mcolor color(arc_x->second, arc_y->second, Mcolor::Intersection);
      if (color.size() > 0) { 
	equiv.addrel(std::make_pair(*arc_x, 0), std::make_pair(*arc_y, 1));
      } 
    } 
  }
  
  equiv.update();
  std::map<colacr_t, std::set<colacr_t> > classes = equiv.get_eclasses<std::set<colacr_t> >(); 

  /*for (const auto &color_set : classes) { 
    std::cerr << color_set.first.first.first << " start new color " << genome_match::mcolor_to_name(color_set.first.first.second) << " in " << color_set.first.second  << std::endl; 
    for (const auto &color : color_set.second) { 
      if (color != color_set.first) {
	std::cerr << color.first.first << " in " << color.second << " " << genome_match::mcolor_to_name(color.first.second) << std::endl;
      } 
    }
    std::cerr << std::endl;
    }*/
  
  size_t count_U = 0; 
  size_t count_cases = classes.size();  
  for (const auto &color_set : classes) { 
    std::map<vertex_t, Mcolor> left; 
    std::map<vertex_t, Mcolor> right; 
    for (const auto &color : color_set.second) { 
      if (color.second == 0) { 
	left.insert(color.first);
      } else { 
	right.insert(color.first);
      }       
    }
    count_U += (left.size() + right.size() - 2);
  }

  return (count_U + count_cases); 
}  

template<class graph_t>
void Algorithm<graph_t>::find_minimal_matching(size_t size_matching, size_t count_vertex, std::set<arc_t>& current, std::set<vertex_t>& processed, const std::map<arc_t, size_t>& weight_edges, std::set<std::set<arc_t> >& answer) {
  for (const auto &edge : weight_edges) { 
    if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) { 
	processed.insert(edge.first.first);
	processed.insert(edge.first.second);
	current.insert(edge.first);
	find_minimal_matching(size_matching + edge.second, count_vertex, current, processed, weight_edges, answer);
	current.erase(edge.first);
	processed.erase(edge.first.first); 
	processed.erase(edge.first.second);
    }
  } 

  //std::cerr << processed.size() << std::endl;
  if (processed.size() == count_vertex) { 
	//std::cerr << "Insert new matching" << std::endl; 
	answer.insert(current);
  } 
} 

template<class graph_t>
std::set<arc_t> Algorithm<graph_t>::create_minimal_matching(const std::set<vertex_t>& vertex_set) {
  std::map<arc_t, size_t> weight_edges; 

  for(const auto& v : vertex_set) { 
    Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(v); 
    Mularcs<Mcolor> mularcs_x = graph->get_adjacent_multiedges(v, true, false); 
    for (const auto& arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
	Mularcs<Mcolor> mularcs_y = graph->get_adjacent_multiedges(arc.first, true, false);
	mularcs_y.erase(v);
	//std::cerr << "Calculate cost " << v << " " << arc.first << " have " << calculate_cost(arc.first, mularcs_x, mularcs_y) << std::endl;
	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), calculate_cost(arc.first, mularcs_x, mularcs_y)));
      }
    } 
  } 

  std::set<vertex_t> temp_processed;
  std::set<arc_t> start_answer;
  for (const auto &edge : weight_edges) { 
	if (postponed_deletions.count(edge.first) != 0) { 
		start_answer.insert(edge.first); 
		temp_processed.insert(edge.first.first);
		temp_processed.insert(edge.first.second);
	} 
  } 

  std::set<std::set<arc_t> > answer;
  for (const auto &edge : weight_edges) { 
    std::set<arc_t> current = start_answer; 
    std::set<vertex_t> processed = temp_processed;  
    if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) { 	
	current.insert(edge.first); 
	processed.insert(edge.first.first); 
	processed.insert(edge.first.second);
    }
    find_minimal_matching(edge.second, vertex_set.size(), current, processed, weight_edges, answer);
  }
 
  std::vector<std::pair<std::set<arc_t>, size_t> > matching; 
  size_t min = 100000; //FIXME
  for (const auto &match: answer) { 
    size_t weight = 0; 
    size_t max = 0; 
    for (const auto &edge: match) { 
      weight += weight_edges.find(edge)->second;
      if (max < weight_edges.find(edge)->second) { 
        max = weight_edges.find(edge)->second;
      } 
    } 
    weight -= max;
    if (weight < min) { 
      min = weight; 
      matching.clear(); 
      matching.push_back(std::make_pair(match, weight));
    } else if (weight == min) { 
      matching.push_back(std::make_pair(match, weight));
    } 
  } 
 
  /*for (const auto& match: matching) { 
  	std::cerr << "New matching have total sum " << match.second << std::endl;
	for (const auto &edge : match.first) { 
		std::cerr << edge.first << "-" << edge.second << std::endl;
	} 
  }*/

  return matching[0].first;
} 

template<class graph_t>
size_t Algorithm<graph_t>::process_minimal_matching(const arc_t& matching) { 
  size_t num_rear = 0; 
  
  Mularcs<Mcolor> mularcs_x = graph->get_adjacent_multiedges(matching.first, true, false);
  Mularcs<Mcolor> mularcs_y = graph->get_adjacent_multiedges(matching.second, true, false);

  typedef std::pair<std::pair<vertex_t, Mcolor>, size_t> colacr_t;
  equivalence<colacr_t> equiv; 

  for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
    if (arc_x->first == matching.second) { 
      continue;
    } 
    for (auto arc_y =  mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
      if (arc_y->first == matching.first) {
	continue;
      }
      Mcolor color(arc_x->second, arc_y->second, Mcolor::Intersection);
      if (color.size() > 0) { 
	equiv.addrel(std::make_pair(*arc_x, 0), std::make_pair(*arc_y, 1));
      } 
    } 
  }
  
  //std::cerr << "start work " << matching.first << " " << matching.second << std::endl; 
  equiv.update();
  std::map<colacr_t, std::set<colacr_t> > classes = equiv.get_eclasses<std::set<colacr_t> >(); 
  
  for (const auto &color_set : classes) { 
    std::multimap<vertex_t, Mcolor> left; 
    std::multimap<vertex_t, Mcolor> right; 
    for (const auto &color : color_set.second) { 
      if (color.second == 0) { 
	left.insert(color.first);
      } else { 
	right.insert(color.first);
      }       
    }

    //std::cerr << "go left" << std::endl; 
    for (auto l = (++left.cbegin()); l != left.cend(); ++l) { 
	vertex_t v = graph->get_adjacent_multiedges(left.cbegin()->first, true, false).get_vertex(l->second);
	assert(!v.empty());
	//std::cerr << "Two break " << left.cbegin()->first << " " <<  v << " " << matching.first << " " << l->first << " " << genome_match::mcolor_to_name(l->second) << std::endl;

	TwoBreak<Mcolor> br(left.cbegin()->first, v, matching.first, l->first, l->second);
	graph->apply_two_break(br);
	Mcolor temp(left.cbegin()->second, l->second, Mcolor::Union);
	left.begin()->second = temp;
	++num_rear;
    } 

    //std::cerr << "go right" << std::endl;
    for (auto r = (++right.cbegin()); r != right.cend(); ++r) { 
	vertex_t v = graph->get_adjacent_multiedges(right.cbegin()->first, true, false).get_vertex(r->second);
	assert(!v.empty());
	//std::cerr << "Two break " << right.cbegin()->first << " " <<  v << " " << matching.second << " " << r->first << " " << genome_match::mcolor_to_name(r->second) << std::endl;
	TwoBreak<Mcolor> br(right.cbegin()->first, v, matching.second, r->first, r->second);
	graph->apply_two_break(br);
	Mcolor temp(right.cbegin()->second, r->second, Mcolor::Union);
	right.begin()->second = temp;
	++num_rear;
    } 

    assert(right.cbegin()->second == left.cbegin()->second);
    assert(graph->is_vec_T_consistent_color(right.cbegin()->second));
    TwoBreak<Mcolor> br(matching.first, left.cbegin()->first, matching.second, right.cbegin()->first, left.cbegin()->second);
    graph->apply_two_break(br);
    ++num_rear;
  } 

  return num_rear; 
} 

template<class graph_t>
bool Algorithm<graph_t>::stage6() { 
  size_t number_break = 0;
  equivalence<vertex_t> CC; // connected components

  for(const auto &x : *graph) {
    Mularcs<Mcolor> mularcs = graph->get_adjacent_multiedges(x); 

    if (mularcs.size() == 1 && mularcs.cbegin()->second == graph->get_complete_color()) { 
      continue; // ignore complete multiedges
    } 

    for(auto im = mularcs.cbegin(); im != mularcs.cend(); ++im) {    
      if (im->first != Infty) { 
	CC.addrel(x, im->first);
      } 
    }
  }
		    
  CC.update();
   
  std::map<vertex_t, std::set<vertex_t> > classes = CC.get_eclasses<std::set<vertex_t> >(); 
 
  for (const auto &vertex_set : classes) { 
    /*std::cerr << "Start new component " << vertex_set.second.size() << std::endl; 
    for (const auto &v : vertex_set.second) { 
      std::cerr << v << " - ";
    }
    std::cerr << std::endl;*/
    std::set<arc_t> matching = create_minimal_matching(vertex_set.second);
    number_break += process_minimal_matching(*matching.begin());
  }
  
  return (number_break != 0);
} 

#endif
