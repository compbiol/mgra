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
  std::map<arc_t, std::pair<size_t, mcolor_t> > weight_edges; 

  for(const auto& v : vertex_set) { 
    mularcs_t&& mularcs = graph->get_adjacent_multiedges(v);
    for (const auto& arc: mularcs) {
      if (arc.first != Infty && weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
        const mularcs_t& mularcs_x = graph->get_adjacent_multiedges_with_info(v, false); 
	mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(arc.first, false);
	mularcs_y.erase(v);
	//std::cerr << "Calculate cost " << v << " " << arc.first << " have " << calculate_cost(arc.first, mularcs_x, mularcs_y) << std::endl;
	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(calculate_cost(arc.first, mularcs_x, mularcs_y), arc.second)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<arc_t> current; 
  std::cerr << vertex_set.size() << std::endl;
  for (const auto &edge : weight_edges) { 
    if (postponed_deletions.defined(edge.first.first, edge.first.second) || (graph->is_T_consistent_color(edge.second.second) && !graph->is_vec_T_consistent_color(edge.second.second)) || edge.second.second.includes(graph->get_root_color())) { 
      current.insert(edge.first); 
      processed.insert({edge.first.first, edge.first.second});
    } 
  } 

  std::set<std::set<arc_t> > answer;
  if (processed.size() == vertex_set.size()) { 
    answer.insert(current);
  } else { 
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
  } 

  std::set<arc_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (const auto &match: answer) { 
    size_t weight = 0; 
    auto max_edge = *match.cbegin(); 
    for (const auto &edge: match) { 
      weight += weight_edges.find(edge)->second.first;
      if (weight_edges.find(max_edge)->second.first < weight_edges.find(edge)->second.first) { 
        max_edge = edge;
      } 
    } 
    if ((weight - weight_edges.find(max_edge)->second.first) < min) { 
      min = weight; 
      minimal_matching = match;
      minimal_matching.erase(max_edge);
    } 
  } 
 
  return minimal_matching;
} 

template<class graph_t>
size_t Algorithm<graph_t>::process_minimal_matching(const std::set<arc_t>& matchings) { 
  size_t num_rear = 0; 
  
  for (const auto& matching : matchings) { 
    mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(matching.first, false);
    mularcs_x.erase(matching.second);
    mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(matching.second, false);
    mularcs_y.erase(matching.first);

    std::cerr << "Start process " << matching.first << " " << matching.second << std::endl;
 
    typedef std::pair<std::pair<vertex_t, mcolor_t>, size_t> colacr_t;
    utility::equivalence<colacr_t> equiv; 

    std::cerr << "Element x " << std::endl;
    for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
      std::cerr << "New edge " << arc_x->first << " " << genome_match::mcolor_to_name(arc_x->second) << std::endl;
    } 

    std::cerr << "Element y " << std::endl;
    for (auto arc_y =  mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
      std::cerr << "New edge " << arc_y->first << " " << genome_match::mcolor_to_name(arc_y->second) << std::endl;
    } 

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
    std::cerr << classes.size() << std::endl;
    
    std::vector<twobreak_t> history;  
    bool good = true;  
    for (auto color_set = classes.cbegin(); (color_set != classes.cend()) && good; ++color_set) { 
      std::multimap<vertex_t, mcolor_t> left; 
      std::multimap<vertex_t, mcolor_t> right; 
      for (const auto &color : color_set->second) { 
        if (color.second == 0) {
          left.insert(color.first);
        } else {
          right.insert(color.first);
        } 
      }
	
      std::cerr << left.size() << " " << right.size() << std::endl;
      for (auto l = (++left.cbegin()); (l != left.cend()) && good; ++l) { 
        if (left.cbegin()->first != l->first) {
	  const vertex_t& v = graph->get_adjacent_multiedges_with_info(left.cbegin()->first, false).get_vertex(l->second); 
          if (!v.empty()) { 
      	    history.push_back(twobreak_t(left.cbegin()->first, v, matching.first, l->first, l->second));
          } else { 
            good = false; 
          } 
	} 
	left.begin()->second = mcolor_t(left.cbegin()->second, l->second, mcolor_t::Union); 
      } 

      for (auto r = (++right.cbegin()); (r != right.cend()) && good; ++r) { 
        if (right.cbegin()->first != r->first) {
	  vertex_t v = graph->get_adjacent_multiedges_with_info(right.cbegin()->first, false).get_vertex(r->second);
          if (!v.empty()) { 
            history.push_back(twobreak_t(right.cbegin()->first, v, matching.second, r->first, r->second));
          } else { 
            good = false; 
          } 
        } 
        right.begin()->second = mcolor_t(right.cbegin()->second, r->second, mcolor_t::Union);
      } 
      
      if (good) {
        assert(right.cbegin()->second == left.cbegin()->second);
        assert(graph->is_vec_T_consistent_color(right.cbegin()->second));
        history.push_back(twobreak_t(matching.first, left.cbegin()->first, matching.second, right.cbegin()->first, left.cbegin()->second));    
      } 
    }

    if (good) {
      for (const auto& break2 : history) {
        graph->apply_two_break(break2);
        ++num_rear; 
      }  
    } 

    if (num_rear != 0) {
      break;
    } 
  } 
  return num_rear; 
} 

template<class graph_t>
bool Algorithm<graph_t>::stage6() { 
  size_t number_break = 0;
  const std::map<vertex_t, std::set<vertex_t> >& classes = graph->split_on_components();
  for (const auto &vertex_set : classes) { 
    std::cerr << vertex_set.second.size() << std::endl;
    if (vertex_set.second.size() <= max_size_component) { 
      std::cerr << "start process component" << std::endl;
      for (auto it = vertex_set.second.cbegin(); it != vertex_set.second.cend(); ++it) {
	std::cerr << *it << " "; 
      } 
      std::cerr << std::endl;
      std::set<arc_t> matching = create_minimal_matching(vertex_set.second);
      std::cerr << "find minimal matching " << matching.size() << std::endl;
      for (auto it = matching.cbegin(); it != matching.cend(); ++it) {
	std::cerr << "(" << it->first  << ", " << it->second << ")" << std::endl; 
      } 
      std::cerr << std::endl;
  
      number_break += process_minimal_matching(matching);
    } 
  }  
  return (number_break != 0);
} 

#endif
