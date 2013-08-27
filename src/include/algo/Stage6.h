#ifndef STAGE6_H_ 
#define STAGE6_H_ 

template<class graph_t>
size_t Algorithm<graph_t>::calculate_cost(const vertex_t& y, const mularcs_t& mularcs_x, const mularcs_t& mularcs_y) { 
  typedef std::pair<std::pair<vertex_t, mcolor_t>, size_t> colacr_t;
  utility::equivalence<colacr_t> equiv; 

  for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
    if (arc_x->first != y) { 
      for (auto arc_y = mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
        mcolor_t color(arc_x->second, arc_y->second, mcolor_t::Intersection);
        if (color.size() > 0) { 
	  equiv.addrel(std::make_pair(*arc_x, 0), std::make_pair(*arc_y, 1));
        } 
      }
    } 
  }
  
  equiv.update();
  const std::map<colacr_t, std::set<colacr_t> >& classes = equiv.get_eclasses<std::set<colacr_t> >(); 

  size_t count_U = 0; 
  for (const auto &color_set : classes) { 
    std::unordered_set<vertex_t> left; 
    std::unordered_set<vertex_t> right;
    for (const auto &color : color_set.second) { 
      (color.second == 0)?left.insert(color.first.first):right.insert(color.first.first);
    }
    if (left.size() != 0) { 
      count_U += (left.size() - 1);
    }  
    if (right.size() != 0) { 
      count_U += (right.size() - 1);
    } 
  }

  return (count_U + classes.size()); 
}  

template<class graph_t>
std::set<arc_t> Algorithm<graph_t>::create_minimal_matching(const std::set<vertex_t>& vertex_set) {
  std::map<arc_t, size_t> weight_edges; 

  for(const auto& v : vertex_set) { 
    mularcs_t mularcs = graph->get_adjacent_multiedges(v); //FIXME
    const mularcs_t& mularcs_x = graph->get_adjacent_multiedges_with_info(v, true, false, false); 
    for (const auto& arc: mularcs) {
      if (arc.first != Infty && weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
	mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(arc.first, true, false, false);
	mularcs_y.erase(v);
	//std::cerr << "Calculate cost " << v << " " << arc.first << " have " << calculate_cost(arc.first, mularcs_x, mularcs_y) << std::endl;
	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), calculate_cost(arc.first, mularcs_x, mularcs_y)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<arc_t> current;
  for (const auto &edge : weight_edges) { 
    if (postponed_deletions.count(edge.first) != 0) { 
      current.insert(edge.first); 
      processed.insert({edge.first.first, edge.first.second});
    } 
  } 

  std::set<std::set<arc_t> > answer;
  std::function<void()> findLambda = [&] () -> void { 
    for (const auto &e : weight_edges) { 
      if (processed.count(e.first.first) == 0 && processed.count(e.first.second) == 0) { 
        processed.insert({e.first.first, e.first.second});
	current.insert(e.first);
	findLambda();
	current.erase(e.first);
	processed.erase(e.first.first);
	processed.erase(e.first.second);
      }
    }
    if (processed.size() == vertex_set.size()) { 
      answer.insert(current);
    } 
  }; 
 
  for (const auto &edge : weight_edges) { 
    if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) { 	
      current.insert(edge.first); 
      processed.insert({edge.first.first, edge.first.second});
      findLambda();  
      current.erase(edge.first);
      processed.erase(edge.first.first);
      processed.erase(edge.first.second);
    }
  }
 
  std::set<arc_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (const auto &match: answer) { 
    size_t weight = 0; 
    auto max_edge = *match.cbegin(); 
    for (const auto &edge: match) { 
      weight += weight_edges.find(edge)->second;
      if (weight_edges.find(max_edge)->second < weight_edges.find(edge)->second) { 
        max_edge = edge;
      } 
    } 
    if ((weight - weight_edges.find(max_edge)->second) < min) { 
      min = weight; 
      minimal_matching = match;
      minimal_matching.erase(max_edge);
    } 
  } 
 
  return minimal_matching;
} 

template<class graph_t>
size_t Algorithm<graph_t>::process_minimal_matching(const arc_t& matching) { 
  size_t num_rear = 0; 
  
  mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(matching.first, true, false, false);
  mularcs_x.erase(matching.second);
  mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(matching.second, true, false, false);
  mularcs_y.erase(matching.first);

  //std::cerr << "Start process " << matching.first << " " << matching.second << std::endl;

  typedef std::pair<std::pair<vertex_t, mcolor_t>, size_t> colacr_t;
  utility::equivalence<colacr_t> equiv; 

  for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
    for (auto arc_y =  mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
      mcolor_t color(arc_x->second, arc_y->second, mcolor_t::Intersection);
      if (color.size() > 0) { 
	equiv.addrel(std::make_pair(*arc_x, 0), std::make_pair(*arc_y, 1));
      } 
    } 
  }
  
  equiv.update();
  const std::map<colacr_t, std::set<colacr_t> >& classes = equiv.get_eclasses<std::set<colacr_t> >(); 
  
  for (const auto &color_set : classes) { 
    std::multimap<vertex_t, mcolor_t> left; 
    std::multimap<vertex_t, mcolor_t> right; 
    for (const auto &color : color_set.second) { 
      if (color.second == 0) {
        left.insert(color.first);
      } else {
        right.insert(color.first);
      } 
    }
	
    //std::cerr << "go left" << std::endl;
    for (auto l = (++left.cbegin()); l != left.cend(); ++l) { 
        if (left.cbegin()->first != l->first) {
	  vertex_t v = graph->get_adjacent_multiedges_with_info(left.cbegin()->first, true, false, false).get_vertex(l->second); 
          assert(!v.empty());
	  graph->apply_two_break(twobreak_t(left.cbegin()->first, v, matching.first, l->first, l->second));
          ++num_rear;
        } 
	left.begin()->second = mcolor_t(left.cbegin()->second, l->second, mcolor_t::Union); 
	
    } 

    //std::cerr << "go right" << std::endl;
    for (auto r = (++right.cbegin()); r != right.cend(); ++r) { 
      //std::cerr << right.cbegin()->first << " " << r->first << " " << genome_match::mcolor_to_name(r->second) << std::endl;
      if (right.cbegin()->first != r->first) {
	vertex_t v = graph->get_adjacent_multiedges_with_info(right.cbegin()->first, true, false, false).get_vertex(r->second);
	assert(!v.empty()); 
        //std::cerr << right.cbegin()->first << " " << v << " " << matching.second << " " << r->first << " " << genome_match::mcolor_to_name(r->second) << std:: endl;
	graph->apply_two_break(twobreak_t(right.cbegin()->first, v, matching.second, r->first, r->second));
	++num_rear;
      } 
      right.begin()->second = mcolor_t(right.cbegin()->second, r->second, mcolor_t::Union);
    } 
    //std::cerr << "finally" << std::endl;
    assert(right.cbegin()->second == left.cbegin()->second);
    assert(graph->is_vec_T_consistent_color(right.cbegin()->second));
    graph->apply_two_break(twobreak_t(matching.first, left.cbegin()->first, matching.second, right.cbegin()->first, left.cbegin()->second));
    ++num_rear;
  } 

  return num_rear; 
} 

template<class graph_t>
bool Algorithm<graph_t>::stage6() { 
  size_t number_break = 0;
  const std::map<vertex_t, std::set<vertex_t> >& classes = graph->split_on_components();
  for (const auto &vertex_set : classes) { 
    if (vertex_set.second.size() <= max_size_component) { 
      std::set<arc_t> matching = create_minimal_matching(vertex_set.second);
      number_break += process_minimal_matching(*matching.begin());
    } 
  }  
  return (number_break != 0);
} 

#endif
