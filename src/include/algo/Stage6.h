#ifndef STAGE6_H_ 
#define STAGE6_H_ 

template<class graph_t>
size_t Algorithm<graph_t>::calculate_cost(vertex_t const & y, mularcs_t const & mularcs_x, mularcs_t const & mularcs_y) { 
  if (y == Infty) {
    return mularcs_x.size(); 
  } 

  typedef std::tuple<vertex_t, mcolor_t, size_t> colacr_t; 
  utility::equivalence<colacr_t> equiv; 
  
  for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
    if (arc_x->first != y) { 
      for (auto arc_y = mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
        mcolor_t color(arc_x->second, arc_y->second, mcolor_t::Intersection);
        if (color.size() > 0) { 
	  equiv.addrel(std::make_tuple(arc_x->first, arc_x->second, 0), std::make_tuple(arc_y->first, arc_y->second, 1));
        } 
      }
    } 
  }
  
  equiv.update();
  std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.template get_eclasses<std::set<colacr_t> >(); 

  size_t count_U = 0; 
  for (auto const & color_set : classes) { 
    count_U += (color_set.second.size() - 1);
  }

  return count_U;
}  

template<class graph_t>
std::multimap<size_t, arc_t> Algorithm<graph_t>::create_minimal_matching(std::set<vertex_t> const & vertex_set) {
  std::map<arc_t, std::pair<size_t, mcolor_t> > weight_edges; 

  for(auto const & v : vertex_set) { 
    mularcs_t mularcs = graph->get_adjacent_multiedges(v);
    for (auto const & arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
        mularcs_t const & mularcs_x = graph->get_adjacent_multiedges_with_info(v, false); 
	
	mularcs_t mularcs_y;
        if (arc.first != Infty) {
          mularcs_y = graph->get_adjacent_multiedges_with_info(arc.first, false);
	  mularcs_y.erase(v);
        } 
        size_t sz = calculate_cost(arc.first, mularcs_x, mularcs_y);
	//std::cerr << "Calculate cost " << v << " " << arc.first << " have " << sz << std::endl;
	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(sz, arc.second)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<arc_t> current; 
  //std::cerr << "Size " << vertex_set.size() << std::endl;
  for (auto const &edge : weight_edges) { 
    if (postponed_deletions.defined(edge.first.first, edge.first.second) || (graph->is_T_consistent_color(edge.second.second) && !graph->is_vec_T_consistent_color(edge.second.second)) || edge.second.second.includes(graph->get_root_color())) { 
      current.insert(edge.first); 
      processed.insert({edge.first.first, edge.first.second});
      processed.erase(Infty);
    } 
  } 

  std::set<std::set<arc_t> > answer;
  if (processed.size() == vertex_set.size()) { 
    answer.insert(current);
  } else { 
    std::function<void()> findLambda = [&] () -> void { 
      for (auto const &e : weight_edges) { 
        if (processed.count(e.first.first) == 0 && processed.count(e.first.second) == 0) { 
          processed.insert({e.first.first, e.first.second});
          processed.erase(Infty);
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
 
    for (auto const &edge : weight_edges) { 
      if (processed.count(edge.first.first) == 0 && processed.count(edge.first.second) == 0) { 	
        current.insert(edge.first); 
        processed.insert({edge.first.first, edge.first.second});
        processed.erase(Infty);
        findLambda();  
        current.erase(edge.first);
        processed.erase(edge.first.first);
        processed.erase(edge.first.second);
      }
    }
  } 

  std::multimap<size_t, arc_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (auto const & match: answer) { 
    size_t weight = 0; 
    //auto max_edge = *match.cbegin(); 
    for (auto const &edge: match) { 
      weight += weight_edges.find(edge)->second.first;
      //if (weight_edges.find(max_edge)->second.first < weight_edges.find(edge)->second.first) { 
        //max_edge = edge;
      //} 
    } 
    if ((weight/* - weight_edges.find(max_edge)->second.first*/) < min) { 
      min = weight; 
      for (auto const &edge: match) { 
        minimal_matching.insert(std::make_pair(weight_edges.find(edge)->second.first, edge));
      } 
      //minimal_matching.erase(max_edge);
    } 
  } 
 
  /*std::cerr << "Find minimal matching " << minimal_matching.size() << std::endl;
  for (auto const & edge : minimal_matching) {
    std::cerr << "(" << edge.second.first  << ", " << edge.second.second << ") " << edge.first << std::endl; 
  }*/

  return minimal_matching;
} 


template<class graph_t>
size_t Algorithm<graph_t>::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  size_t num_rear = 0;

  //std::cerr << "Start to work with (" << x << ", " << y << ") " << genome_match::mcolor_to_name(color) << std::endl;  

  if (y == Infty) {
    mcolor_t need_color(color, graph->get_adjacent_multiedges(x).get_multicolor(y), mcolor_t::Difference);    
    mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    for (auto const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && graph->is_vec_T_consistent_color(arc.second)) {
	graph->apply_two_break(twobreak_t(x, arc.first, Infty, Infty, arc.second));
	++num_rear; 
      } 
    } 
  } else { 
    mcolor_t need_color(color, graph->get_adjacent_multiedges(x).get_multicolor(y), mcolor_t::Difference);
    //std::cerr << "Target color " << genome_match::mcolor_to_name(need_color) << " " << graph->is_vec_T_consistent_color(need_color) << std::endl; 
    //std::cerr << " edge: " << genome_match::mcolor_to_name(graph->get_adjacent_multiedges(x).get_multicolor(y)) << std::endl;
    mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

    typedef std::pair<std::pair<vertex_t, structure::Mcolor>, size_t> colacr_t;
    utility::equivalence<colacr_t> equiv; 

    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) { 
        mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
        if (color.size() > 0) { 
	  equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
        } 
      } 
    }  
    equiv.update();

    std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.get_eclasses<std::set<colacr_t> >();   
    for (auto const & color_set : classes) { 
      std::multimap<mcolor_t, vertex_t> left; 
      std::multimap<mcolor_t, vertex_t> right; 
      for (auto const & color : color_set.second) { 
        if (color.second == 0) {
          left.insert(std::make_pair(color.first.second, color.first.first));
        } else {
          right.insert(std::make_pair(color.first.second, color.first.first));
        } 
      }
	
      if (left.size() == 1 && right.size() == 1 && left.begin()->second == right.begin()->second) {
        assert(graph->is_vec_T_consistent_color(left.begin()->first)); 
        assert(graph->is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          //std::cerr << "Good 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
	  graph->apply_two_break(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
	  ++num_rear; 
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          //std::cerr << "Left and go recursevly" << std::endl;
          num_rear += take_edge_on_color(y, left.begin()->first, right.begin()->second);
          assert(graph->is_vec_T_consistent_color(left.begin()->first)); 
          //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          graph->apply_two_break(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
          ++num_rear; 
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          //std::cerr << "Right and go recursevly" << std::endl;
          num_rear += take_edge_on_color(x, right.begin()->first, left.begin()->second);
	  assert(graph->is_vec_T_consistent_color(right.begin()->first));  
          //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(right.begin()->first) << std::endl;
          graph->apply_two_break(twobreak_t(x, left.begin()->second, y, right.begin()->second, right.begin()->first));
          ++num_rear; 
        }   
      } else { 
        assert(left.size() == 1 || right.size() == 1); 
      } 
    } 
  } 

  return num_rear;
} 


/*template<class graph_t>
size_t Algorithm<graph_t>::process_minimal_matching(std::set<arc_t> const & matchings) { 
  size_t num_rear = 0; 
  
  for (auto const & matching : matchings) { 
    mularcs_t mularcs_x = graph->get_adjacent_multiedges_with_info(matching.first, false);
    mularcs_x.erase(matching.second);
    mularcs_t mularcs_y = graph->get_adjacent_multiedges_with_info(matching.second, false);
    mularcs_y.erase(matching.first);

    std::cerr << "Start process " << matching.first << " " << matching.second << std::endl;
 
    typedef std::pair<std::pair<vertex_t, structure::Mcolor>, size_t> colacr_t;
    utility::equivalence<colacr_t> equiv; 

    //std::cerr << "Element x " << std::endl;
    //for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
    //  std::cerr << "New edge " << arc_x->first << " " << genome_match::mcolor_to_name(arc_x->second) << std::endl;
    //} 

    //std::cerr << "Element y " << std::endl;
    //for (auto arc_y =  mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
    //  std::cerr << "New edge " << arc_y->first << " " << genome_match::mcolor_to_name(arc_y->second) << std::endl;
    //} 

    for (auto arc_x = mularcs_x.cbegin(); arc_x != mularcs_x.cend(); ++arc_x) { 
      for (auto arc_y =  mularcs_y.cbegin(); arc_y != mularcs_y.cend(); ++arc_y) { 
        mcolor_t color(arc_x->second, arc_y->second, mcolor_t::Intersection);
        if (color.size() > 0) { 
	  equiv.addrel(std::make_pair(*arc_x, 0), std::make_pair(*arc_y, 1));
        } 
      } 
    }
  
    equiv.update();
    std::map<colacr_t, std::set<colacr_t> > const & classes = equiv.get_eclasses<std::set<colacr_t> >();   
    std::cerr << classes.size() << std::endl;
    
    std::vector<twobreak_t> history;  
    bool good = true;  
    for (auto color_set = classes.cbegin(); (color_set != classes.cend()) && good; ++color_set) { 
      std::multimap<vertex_t, mcolor_t> left; 
      std::multimap<vertex_t, mcolor_t> right; 
      for (auto const & color : color_set->second) { 
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
*/

template<class graph_t>
bool Algorithm<graph_t>::stage6() { 
  size_t number_break = 0;

  std::map<vertex_t, std::set<vertex_t> > const & classes = graph->split_on_components();
  for (auto const & vertex_set : classes) { 
 
    std::cerr << "Component have size " << vertex_set.second.size() << std::endl;

    if (vertex_set.second.size() <= max_size_component) { 

      //std::cerr << "Start process component" << std::endl;
      /*for (auto const v : vertex_set.second) {
	std::cerr << v << " "; 
      } 
      std::cerr << std::endl;*/

      std::multimap<size_t, arc_t> matching = create_minimal_matching(vertex_set.second);

      //std::cerr << "Start process minimal matching" << std::endl;
      auto edge = matching.cbegin()->second;
      number_break += take_edge_on_color(edge.first, graph->get_complete_color(), edge.second);
    } 
  }  

  return (number_break != 0);
} 

#endif
