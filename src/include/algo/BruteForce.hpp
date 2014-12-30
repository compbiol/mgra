#ifndef BRUTE_FORCE_HPP 
#define BRUTE_FORCE_HPP

template<class graph_t>
struct Algorithm<graph_t>::BruteForce : public Algorithm<graph_t>::Stage {
  typedef Stage base;

  typedef typename graph_t::mcolor_type mcolor_t;
  typedef typename graph_t::edge_t edge_t;
  typedef typename graph_t::mularcs_t mularcs_t; 
  typedef typename graph_t::twobreak_t twobreak_t;
  
  BruteForce(std::shared_ptr<graph_t> const & graph, size_t size_component)
  : Stage(graph) 
  , max_size_component(size_component)
  {
  }
    
  bool do_action() override;

  std::string get_name() override {
    return "Resolve component with score function."; 
  }

  size_t take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y);
  
private:
  std::multimap<size_t, edge_t> create_minimal_matching(std::set<vertex_t> const & vertex_set); 
  
private: 
  size_t max_size_component;
};

template<class graph_t>
bool Algorithm<graph_t>::BruteForce::do_action() { 
  INFO("Run BruteForce stage")
  size_t number_rear = 0; // number of rearrangements 

  utility::equivalence<vertex_t> components = this->graph->split_on_components();
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();

  for (auto const & vertex_set : classes) { 
     if ((vertex_set.second.size() <= max_size_component)) { 
      std::multimap<size_t, edge_t> matching;
      matching = create_minimal_matching(vertex_set.second);
      auto edge = matching.cbegin()->second;
      number_rear += take_edge_on_color(edge.first, this->graph->get_complete_color(), edge.second);
    } 
  }

  return (number_rear != 0);
} 

template<class graph_t>
std::multimap<size_t, typename graph_t::edge_t> Algorithm<graph_t>::BruteForce::create_minimal_matching(std::set<vertex_t> const & vertex_set) {
  std::map<edge_t, std::pair<size_t, mcolor_t> > weight_edges; 

  for(vertex_t const & v : vertex_set) { 
    mularcs_t mularcs = this->graph->get_all_adjacent_multiedges(v);
    for (auto const & arc: mularcs) {
      if (weight_edges.count(std::make_pair(v, arc.first)) == 0 && weight_edges.count(std::make_pair(arc.first, v)) == 0) { 
        size_t sz = this->graph->calculate_cost(v, arc.first);
      	//std::cerr << "Calculate cost " << v << " " << arc.first << " have " << sz << std::endl;
      	weight_edges.insert(std::make_pair(std::make_pair(v, arc.first), std::make_pair(sz, arc.second)));
      }
    } 
  } 

  std::unordered_set<vertex_t> processed;
  std::set<edge_t> current; 
  for (auto const &edge : weight_edges) { 
    bool is_have_tc = false; 
    auto edge_colors = this->graph->split_color(edge.second.second);
    
    for (auto const & color : edge_colors) { 
      if (this->graph->is_T_consistent_color(color) && !this->graph->is_vec_T_consistent_color(color)) { 
        is_have_tc = true; 
        break;
      }
    }

    if (this->graph->is_postponed_deletion(edge.first.first, edge.first.second) || is_have_tc) { 
      current.insert(edge.first); 
      processed.insert(edge.first.first); processed.insert(edge.first.second);
      processed.erase(Infty);
    } 
  } 

  std::set<std::set<edge_t> > answer;
  if (processed.size() == vertex_set.size()) { 
    answer.insert(current);
  } else { 
    std::function<void()> findLambda = [&] () -> void { 
      for (auto const &e : weight_edges) { 
        if (processed.count(e.first.first) == 0 && processed.count(e.first.second) == 0) { 
          processed.insert(e.first.first); processed.insert(e.first.second);
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
        processed.insert(edge.first.first); processed.insert(edge.first.second);
        processed.erase(Infty);
        findLambda();  
        current.erase(edge.first);
        processed.erase(edge.first.first);
        processed.erase(edge.first.second);
      }
    }
  } 

  std::multimap<size_t, edge_t> minimal_matching; 
  size_t min = std::numeric_limits<size_t>::max() / 4;
  for (auto const & match: answer) { 
    size_t weight = 0; 
    for (auto const &edge: match) { 
      weight += weight_edges.find(edge)->second.first;
    } 
    if (weight < min) { 
      min = weight; 
      minimal_matching.clear();
      for (auto const &edge: match) { 
        minimal_matching.insert(std::make_pair(weight_edges.find(edge)->second.first, edge));
      } 
    } 
  } 
 
  /*std::cerr << "Find minimal matching " << minimal_matching.size() << std::endl;
  for (auto const & edge : minimal_matching) {
    std::cerr << "(" << edge.second.first  << ", " << edge.second.second << ") " << edge.first << std::endl; 
  }*/

  return minimal_matching;
} 


template<class graph_t>
size_t Algorithm<graph_t>::BruteForce::take_edge_on_color(vertex_t const & x, mcolor_t const & color, vertex_t const & y) {
  size_t num_rear = 0;

  //std::cerr << "Start to work with (" << x << ", " << y << ") " << genome_match::mcolor_to_name(color) << std::endl;  
  std::list<twobreak_t> result; 

  if (y == Infty) {
    mcolor_t need_color(color, this->graph->get_all_multicolor_edge(x, y), mcolor_t::Difference);    
    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    for (auto const & arc : mularcs_x) {
      if (need_color.includes(arc.second) && this->graph->is_vec_T_consistent_color(arc.second)) {
        twobreak_t twobreak(x, arc.first, Infty, Infty, arc.second);
      	this->graph->apply(twobreak);
        result.push_back(twobreak);
      	++num_rear; 
      } 
    } 
  } else { 
    mcolor_t need_color(color, this->graph->get_all_multicolor_edge(x, y), mcolor_t::Difference);
    //std::cerr << "Target color " << genome_match::mcolor_to_name(need_color) << " " << this->graph->is_vec_T_consistent_color(need_color) << std::endl; 
    //std::cerr << " edge: " << genome_match::mcolor_to_name(this->graph->get_all_multicolor_edge(x, y)) << std::endl;
  
    mularcs_t mularcs_x = this->graph->get_all_adjacent_multiedges_with_info(x, false);
    mularcs_x.erase(y);
  
    mularcs_t mularcs_y = this->graph->get_all_adjacent_multiedges_with_info(y, false);
    mularcs_y.erase(x);

    /*std::cerr << "mulacrs_x vertex have " << std::endl; 
    for (auto const & l : mularcs_x) { 
      std::cerr << genome_match::mcolor_to_name(l.second) << " " << this->graph->is_vec_T_consistent_color(l.second) << " " << l.first << std::endl;
    }

    std::cerr << "mularcs_y vertex have " << std::endl; 
    for (auto const & r : mularcs_y) {
      std::cerr << genome_match::mcolor_to_name(r.second) << " " << this->graph->is_vec_T_consistent_color(r.second) << " " << r.first << std::endl;
    }*/
  

    typedef std::pair<std::pair<vertex_t, structure::Mcolor>, size_t> colacr_t;
    utility::equivalence<colacr_t> equiv; 

    for (auto const & arc_x : mularcs_x) { 
      for (auto const & arc_y : mularcs_y) {
        if (need_color.includes(arc_x.second) && need_color.includes(arc_y.second)) {  
          mcolor_t color(arc_x.second, arc_y.second, mcolor_t::Intersection);
          if (color.size() > 0) {   
            equiv.addrel(std::make_pair(arc_x, 0), std::make_pair(arc_y, 1));
          }
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
	
      
      /*
      std::cerr << "left vertex have " << std::endl; 
      for (auto const & l : left) { 
        std::cerr << genome_match::mcolor_to_name(l.first) << " " << this->graph->is_vec_T_consistent_color(l.first) << " " << l.second << std::endl;
      }

      std::cerr << "right vertex have " << std::endl; 
      for (auto const & r : right) {
        std::cerr << genome_match::mcolor_to_name(r.first) << " " << this->graph->is_vec_T_consistent_color(r.first) << " " << r.second << std::endl;
      }
      */

      if (left.size() == 1 && right.size() == 1 && left.begin()->first == right.begin()->first) {
        assert(this->graph->is_vec_T_consistent_color(left.begin()->first)); 
        assert(this->graph->is_vec_T_consistent_color(right.begin()->first)); 
        if (need_color.includes(left.begin()->first)) {
          //std::cerr << "Good 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
      	  this->graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
      	  ++num_rear; 
        } 
      } else if (left.size() == 1 || right.size() == 1) {  
        if (left.size() == 1 && need_color.includes(left.begin()->first)) {
          //std::cerr << "Left and go recursevly right" << std::endl;
          num_rear += take_edge_on_color(y, left.begin()->first, right.begin()->second);
          assert(this->graph->is_vec_T_consistent_color(left.begin()->first)); 
          //std::cerr << "Left 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(left.begin()->first) << std::endl;
          this->graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, left.begin()->first));
          ++num_rear; 
        } else if (right.size() == 1 && need_color.includes(right.begin()->first)) {
          //std::cerr << "Right and go recursevly left" << std::endl;
          num_rear += take_edge_on_color(x, right.begin()->first, left.begin()->second);
      	  assert(this->graph->is_vec_T_consistent_color(right.begin()->first));  
          //std::cerr << "Right 2-break: " << x << " " << left.begin()->second << " " << y << " " << right.begin()->second <<  " " << genome_match::mcolor_to_name(right.begin()->first) << std::endl;
          this->graph->apply(twobreak_t(x, left.begin()->second, y, right.begin()->second, right.begin()->first));
          ++num_rear; 
        } else { 
          std::cerr << "Se are here " << std::endl;
        }   
      } else { 
        assert(left.size() == 1 || right.size() == 1); 
      } 
    } 
  } 

  return num_rear;
}

#endif
